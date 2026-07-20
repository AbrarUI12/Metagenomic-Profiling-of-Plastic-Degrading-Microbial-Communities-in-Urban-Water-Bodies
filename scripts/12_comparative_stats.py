"""Comparative statistics across sites and seasons.

Answers the three questions the design was built for, and writes every result as a
CSV so the report and figures never hardcode a number:

1. Does composition differ by **season** or by **site**?  -> PCoA + PERMANOVA
2. Is there an upstream->downstream **gradient**?         -> Spearman vs latitude
3. **Which** families/genera drive the differences?       -> per-feature tests + BH FDR

PCoA (classical MDS) and PERMANOVA are implemented directly on numpy/scipy rather
than pulling in scikit-bio, keeping the environment small and reproducible.

Statistical note: with 18 samples over 8 sites, several sites have a single sample.
The primary two-factor test is therefore run on the **balanced core** (sites sampled
in all three seasons); the full-18 test is reported as secondary.
"""

import argparse
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import kruskal, mannwhitneyu, spearmanr

from scripts.utils import ensure_dir


def bray_curtis(matrix):
    """Square Bray-Curtis distance matrix from a samples x features table."""
    values = np.asarray(matrix, dtype=float)
    # Rows that are entirely zero make Bray-Curtis undefined; nudge them so the
    # distance stays defined (an all-zero sample is maximally distant from any
    # sample that has hits, which is the correct interpretation).
    if values.sum() == 0:
        return np.zeros((values.shape[0], values.shape[0]))
    return squareform(pdist(values, metric="braycurtis"))


def pcoa(distance, n_axes=3):
    """Classical multidimensional scaling (principal coordinates analysis).

    Returns (coordinates, proportion of variance explained per axis).
    """
    n = distance.shape[0]
    a = -0.5 * distance**2
    centering = np.eye(n) - np.ones((n, n)) / n
    gram = centering @ a @ centering
    eigvals, eigvecs = np.linalg.eigh(gram)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    positive = eigvals > 0
    explained = np.zeros_like(eigvals)
    if positive.any():
        explained[positive] = eigvals[positive] / eigvals[positive].sum()
    k = min(n_axes, n - 1)
    coords = eigvecs[:, :k] * np.sqrt(np.clip(eigvals[:k], 0, None))
    return coords, explained[:k]


def permanova(distance, labels, permutations=999, seed=42):
    """One-factor PERMANOVA (Anderson 2001) via pseudo-F and label permutation."""
    labels = np.asarray(labels)
    n = len(labels)
    groups = np.unique(labels)
    a = len(groups)
    if a < 2 or n <= a:
        return {"pseudo_F": np.nan, "R2": np.nan, "p_value": np.nan,
                "n": n, "groups": a, "permutations": 0}

    def ss_within(lab):
        total = 0.0
        for g in groups:
            idx = np.where(lab == g)[0]
            if len(idx) < 2:
                continue
            sub = distance[np.ix_(idx, idx)]
            total += (sub**2).sum() / (2.0 * len(idx))
        return total

    sst = (distance**2).sum() / (2.0 * n)

    def pseudo_f(lab):
        ssw = ss_within(lab)
        ssa = sst - ssw
        if ssw <= 0:
            return np.nan, np.nan
        f = (ssa / (a - 1)) / (ssw / (n - a))
        return f, ssa / sst if sst > 0 else np.nan

    observed_f, r2 = pseudo_f(labels)
    rng = np.random.default_rng(seed)
    permuted = np.empty(permutations)
    shuffled = labels.copy()
    for i in range(permutations):
        rng.shuffle(shuffled)
        permuted[i], _ = pseudo_f(shuffled)
    valid = permuted[~np.isnan(permuted)]
    p = (np.sum(valid >= observed_f) + 1) / (len(valid) + 1) if len(valid) else np.nan
    return {"pseudo_F": observed_f, "R2": r2, "p_value": p,
            "n": n, "groups": a, "permutations": permutations}


def benjamini_hochberg(pvalues):
    """BH-adjusted p-values (same order as input); NaNs pass through."""
    p = np.asarray(pvalues, dtype=float)
    ok = ~np.isnan(p)
    adjusted = np.full_like(p, np.nan)
    if not ok.any():
        return adjusted
    sub = p[ok]
    order = np.argsort(sub)
    ranked = sub[order]
    m = len(ranked)
    q = ranked * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(q, 0, 1)
    adjusted[ok] = out
    return adjusted


def balanced_core(meta):
    """Runs from sites that were sampled in every season - the cleanest design."""
    n_seasons = meta["season"].nunique()
    counts = meta.groupby("site")["season"].nunique()
    sites = counts[counts == n_seasons].index
    return meta[meta["site"].isin(sites)]


def run_permanova_set(matrix, meta, label):
    """PERMANOVA for season / site / urban on one feature matrix."""
    rows = []
    # Fewer than three samples cannot support a permutation test of any factor.
    if len(meta) < 3:
        return rows
    subset = matrix.loc[meta.index]
    distance = bray_curtis(subset)
    for factor in ["season", "site", "urban"]:
        if meta[factor].nunique() < 2:
            continue
        result = permanova(distance, meta[factor].to_numpy())
        result.update({"dataset": label, "factor": factor})
        rows.append(result)
    return rows


def per_feature_tests(matrix, meta):
    """Kruskal-Wallis across seasons and Mann-Whitney urban vs rest, per feature."""
    rows = []
    for feature in matrix.columns:
        values = matrix[feature]
        season_groups = [values[meta["season"] == s].to_numpy() for s in meta["season"].unique()]
        season_p = np.nan
        if len(season_groups) > 1 and all(len(g) > 0 for g in season_groups):
            if len({tuple(g) for g in season_groups}) > 1 and values.nunique() > 1:
                try:
                    season_p = kruskal(*season_groups).pvalue
                except ValueError:
                    season_p = np.nan

        urban_vals = values[meta["urban"] == "1"].to_numpy()
        rural_vals = values[meta["urban"] != "1"].to_numpy()
        urban_p = np.nan
        if len(urban_vals) > 0 and len(rural_vals) > 0 and values.nunique() > 1:
            try:
                urban_p = mannwhitneyu(urban_vals, rural_vals, alternative="two-sided").pvalue
            except ValueError:
                urban_p = np.nan

        rows.append({
            "feature": feature,
            "mean_cpm": float(values.mean()),
            "mean_urban": float(urban_vals.mean()) if len(urban_vals) else np.nan,
            "mean_non_urban": float(rural_vals.mean()) if len(rural_vals) else np.nan,
            "kruskal_season_p": season_p,
            "mannwhitney_urban_p": urban_p,
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df["kruskal_season_q"] = benjamini_hochberg(df["kruskal_season_p"])
        df["mannwhitney_urban_q"] = benjamini_hochberg(df["mannwhitney_urban_p"])
        df = df.sort_values("kruskal_season_p", na_position="last")
    return df


def gradient_tests(summary, matrix):
    """Spearman correlation of potential against latitude (upstream -> downstream)."""
    lat = summary["lat"].astype(float).to_numpy()
    rows = []
    targets = {"total_hits_cpm": summary["total_hits_cpm"].to_numpy()}
    for family in matrix.columns:
        targets[family] = matrix[family].to_numpy()
    for name, values in targets.items():
        if np.nanstd(values) == 0:
            rho, p = np.nan, np.nan
        else:
            rho, p = spearmanr(lat, values)
        rows.append({"feature": name, "spearman_rho_vs_latitude": rho, "p_value": p})
    df = pd.DataFrame(rows)
    df["q_value"] = benjamini_hochberg(df["p_value"])
    return df


def diversity_tests(summary):
    rows = []
    for metric in ["total_hits_cpm", "genus_richness", "genus_shannon", "family_shannon"]:
        values = summary[metric].astype(float)
        season_groups = [values[summary["season"] == s].to_numpy() for s in summary["season"].unique()]
        p_season = np.nan
        if len(season_groups) > 1 and values.nunique() > 1:
            try:
                p_season = kruskal(*season_groups).pvalue
            except ValueError:
                p_season = np.nan
        urban = values[summary["urban"] == "1"].to_numpy()
        rural = values[summary["urban"] != "1"].to_numpy()
        p_urban = np.nan
        if len(urban) and len(rural) and values.nunique() > 1:
            try:
                p_urban = mannwhitneyu(urban, rural, alternative="two-sided").pvalue
            except ValueError:
                p_urban = np.nan
        rows.append({
            "metric": metric,
            "mean_urban": float(urban.mean()) if len(urban) else np.nan,
            "mean_non_urban": float(rural.mean()) if len(rural) else np.nan,
            "kruskal_season_p": p_season,
            "mannwhitney_urban_p": p_urban,
        })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables-dir", default="outputs/comparative/tables")
    parser.add_argument("--out-dir", default="outputs/comparative/tables")
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    tables = Path(args.tables_dir)
    out_dir = Path(args.out_dir)
    ensure_dir(out_dir)

    summary = pd.read_csv(tables / "sample_summary.tsv", sep="\t", index_col="run_accession")
    summary["urban"] = summary["urban"].astype(str)
    enzyme_cpm = pd.read_csv(tables / "sample_enzyme_matrix_cpm.csv", index_col="run_accession")
    genus_cpm = pd.read_csv(tables / "sample_genus_matrix_cpm.csv", index_col="run_accession")

    core = balanced_core(summary)
    print(f"Samples: {len(summary)} total; balanced core: {len(core)} "
          f"({core['site'].nunique()} sites x {core['season'].nunique()} seasons)")

    # --- Ordination -------------------------------------------------------
    for name, matrix in [("enzyme", enzyme_cpm), ("genus", genus_cpm)]:
        distance = bray_curtis(matrix)
        coords, explained = pcoa(distance)
        coord_df = pd.DataFrame(
            coords,
            index=matrix.index,
            columns=[f"PCo{i+1}" for i in range(coords.shape[1])],
        )
        for col in ["site", "season", "lat", "urban"]:
            coord_df[col] = summary.loc[coord_df.index, col]
        coord_df.to_csv(out_dir / f"pcoa_{name}_coords.csv")
        pd.DataFrame({
            "axis": [f"PCo{i+1}" for i in range(len(explained))],
            "proportion_explained": explained,
        }).to_csv(out_dir / f"pcoa_{name}_explained.csv", index=False)
        pd.DataFrame(distance, index=matrix.index, columns=matrix.index).to_csv(
            out_dir / f"braycurtis_{name}.csv"
        )
        shown = ", ".join(
            f"PCo{i+1}={explained[i]:.1%}" for i in range(min(2, len(explained)))
        )
        print(f"  PCoA ({name}): {shown}")

    # --- PERMANOVA --------------------------------------------------------
    permanova_rows = []
    for name, matrix in [("enzyme", enzyme_cpm), ("genus", genus_cpm)]:
        permanova_rows += run_permanova_set(matrix, core, f"{name}_balanced_core")
        permanova_rows += run_permanova_set(matrix, summary, f"{name}_all_samples")
    permanova_df = pd.DataFrame(permanova_rows)
    permanova_df.to_csv(out_dir / "permanova_results.csv", index=False)
    print("\nPERMANOVA:")
    for _, row in permanova_df.iterrows():
        print(f"  {row['dataset']:26} ~{row['factor']:7} "
              f"F={row['pseudo_F']:.3f} R2={row['R2']:.3f} p={row['p_value']:.4f}")

    # --- Per-feature, gradient, diversity ---------------------------------
    per_feature_tests(enzyme_cpm, summary).to_csv(
        out_dir / "enzyme_family_tests.csv", index=False)
    top_genera = genus_cpm.sum().sort_values(ascending=False).head(50).index
    per_feature_tests(genus_cpm[top_genera], summary).to_csv(
        out_dir / "genus_tests_top50.csv", index=False)
    gradient_tests(summary, enzyme_cpm).to_csv(
        out_dir / "gradient_latitude.csv", index=False)
    diversity_tests(summary).to_csv(
        out_dir / "diversity_tests.csv", index=False)

    print(f"\nStats written to {out_dir}/")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[comparative_stats] Error: {exc}", file=sys.stderr)
        sys.exit(1)
