"""Assemble the comparative results into a markdown report.

Every number is read from the CSVs produced upstream; nothing is hardcoded, so
re-running the pipeline with different parameters yields a correct report without
editing this file. Wording stays conservative: we report differences in *putative
genetic potential* inferred by homology, never activity.
"""

import argparse
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

import pandas as pd

from scripts.utils import ensure_dir, utc_now

SEASON_LABELS = {
    "post_monsoon": "post-monsoon (Sep 2021)",
    "winter": "winter / dry (Dec 2021-Jan 2022)",
    "pre_monsoon": "pre-monsoon (Apr 2022)",
}


def significance(p):
    if pd.isna(p):
        return "not testable"
    if p < 0.001:
        return "p < 0.001"
    return f"p = {p:.3f}"


def design_section(summary):
    lines = ["## Study design", ""]
    lines.append(f"- Samples analysed: **{len(summary)}** shotgun sediment metagenomes")
    lines.append(f"- Sites: **{summary['site'].nunique()}** along the Brahmaputra")
    lines.append(f"- Seasons: **{summary['season'].nunique()}**")
    lines.append(f"- Reads searched per sample: **{int(summary['reads_searched'].iloc[0]):,}** "
                 "(equal-depth subsampling, so hit counts are directly comparable)")
    lines.append("")
    lines.append("| Site | Latitude | Urban | Seasons sampled |")
    lines.append("|---|---|---|---|")
    ordered = summary.groupby("site").agg(
        lat=("lat", "first"), urban=("urban", "first"), n=("season", "nunique")
    ).sort_values("lat", ascending=False)
    for site, row in ordered.iterrows():
        urban = "yes" if str(row["urban"]) == "1" else "no"
        lines.append(f"| {site} | {row['lat']} | {urban} | {int(row['n'])} |")
    lines.append("")
    return lines


def overview_section(summary, enzyme_raw):
    total = int(summary["total_hits"].sum())
    lines = ["## Overall putative plastic-degrading potential", ""]
    lines.append(f"- Total filtered enzyme hits across all samples: **{total:,}**")
    family_totals = enzyme_raw.sum().sort_values(ascending=False)
    lines.append("- Hits by enzyme family:")
    for family, count in family_totals.items():
        share = count / total * 100 if total else 0
        lines.append(f"  - {family}: {int(count):,} ({share:.1f}%)")
    hottest = summary["total_hits_cpm"].idxmax()
    coolest = summary["total_hits_cpm"].idxmin()
    lines.append(
        f"- Highest potential: **{summary.loc[hottest, 'site']}** "
        f"({SEASON_LABELS.get(summary.loc[hottest, 'season'], summary.loc[hottest, 'season'])}), "
        f"{summary.loc[hottest, 'total_hits_cpm']:.2f} CPM"
    )
    lines.append(
        f"- Lowest potential: **{summary.loc[coolest, 'site']}** "
        f"({SEASON_LABELS.get(summary.loc[coolest, 'season'], summary.loc[coolest, 'season'])}), "
        f"{summary.loc[coolest, 'total_hits_cpm']:.2f} CPM"
    )
    lines.append("")
    return lines


def permanova_section(permanova):
    lines = ["## Does composition differ by site or season? (PERMANOVA)", ""]
    lines.append("Bray-Curtis distances on CPM profiles; pseudo-F with 999 label permutations. "
                 "The balanced core (sites sampled in every season) is the primary test; the "
                 "all-sample test is secondary because several sites have a single sample.")
    lines.append("")
    lines.append("| Dataset | Factor | pseudo-F | R2 | p |")
    lines.append("|---|---|---|---|---|")
    for _, row in permanova.iterrows():
        f = "n/a" if pd.isna(row["pseudo_F"]) else f"{row['pseudo_F']:.3f}"
        r2 = "n/a" if pd.isna(row["R2"]) else f"{row['R2']:.3f}"
        p = "n/a" if pd.isna(row["p_value"]) else f"{row['p_value']:.4f}"
        lines.append(f"| {row['dataset']} | {row['factor']} | {f} | {r2} | {p} |")
    lines.append("")
    return lines


def gradient_section(gradient):
    lines = ["## Upstream -> downstream gradient", ""]
    lines.append("Spearman correlation against latitude. Latitude decreases downstream, so a "
                 "**negative** rho means potential increases toward the downstream reaches.")
    lines.append("")
    lines.append("| Feature | Spearman rho | p | q (BH) |")
    lines.append("|---|---|---|---|")
    for _, row in gradient.iterrows():
        rho = "n/a" if pd.isna(row["spearman_rho_vs_latitude"]) else f"{row['spearman_rho_vs_latitude']:.3f}"
        p = "n/a" if pd.isna(row["p_value"]) else f"{row['p_value']:.4f}"
        q = "n/a" if pd.isna(row["q_value"]) else f"{row['q_value']:.4f}"
        lines.append(f"| {row['feature']} | {rho} | {p} | {q} |")
    lines.append("")
    return lines


def feature_section(enzyme_tests, diversity):
    lines = ["## Which families and diversity metrics differ?", ""]
    lines.append("Kruskal-Wallis across seasons and Mann-Whitney for urban (Guwahati/Uzanbazar) "
                 "versus non-urban sites, Benjamini-Hochberg corrected.")
    lines.append("")
    lines.append("| Enzyme family | mean CPM | urban | non-urban | season p (q) | urban p (q) |")
    lines.append("|---|---|---|---|---|---|")
    for _, row in enzyme_tests.iterrows():
        lines.append(
            f"| {row['feature']} | {row['mean_cpm']:.3f} | {row['mean_urban']:.3f} | "
            f"{row['mean_non_urban']:.3f} | {significance(row['kruskal_season_p'])} "
            f"({row['kruskal_season_q']:.3f}) | {significance(row['mannwhitney_urban_p'])} "
            f"({row['mannwhitney_urban_q']:.3f}) |"
        )
    lines.append("")
    lines.append("| Diversity metric | urban | non-urban | season p | urban p |")
    lines.append("|---|---|---|---|---|")
    for _, row in diversity.iterrows():
        lines.append(
            f"| {row['metric']} | {row['mean_urban']:.3f} | {row['mean_non_urban']:.3f} | "
            f"{significance(row['kruskal_season_p'])} | {significance(row['mannwhitney_urban_p'])} |"
        )
    lines.append("")
    return lines


def limitations_section():
    return [
        "## Limitations (thesis-safe)",
        "",
        "- All results are **putative genetic potential** inferred by homology. Presence of a "
        "gene fragment is not expression, and not evidence of active plastic degradation.",
        "- Best-hit taxonomy can misassign conserved hydrolases across taxa; genus labels are "
        "indicative, not definitive.",
        "- The reference database is curated from public annotations and spans broad enzyme "
        "families; functional specificity requires experimental validation.",
        "- Bulk sediment metagenomes do not demonstrate plastisphere association.",
        "- PERMANOVA tests differences in composition; it does not establish ecological "
        "interaction, causation, or a mechanism.",
        "- Several sites contribute a single sample, so site-level effects outside the balanced "
        "core are descriptive rather than inferential.",
        "- Short reads (~150 bp) limit resolution; genes are detected as fragments, not "
        "assembled full-length sequences.",
        "",
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables-dir", default="outputs/comparative/tables")
    parser.add_argument("--fig-dir", default="outputs/comparative/figures")
    parser.add_argument("--out-md", default="outputs/comparative/report.md")
    args = parser.parse_args()

    tables = Path(args.tables_dir)
    summary = pd.read_csv(tables / "sample_summary.tsv", sep="\t", index_col="run_accession")
    enzyme_raw = pd.read_csv(tables / "sample_enzyme_matrix_raw.csv", index_col="run_accession")
    permanova = pd.read_csv(tables / "permanova_results.csv")
    gradient = pd.read_csv(tables / "gradient_latitude.csv")
    enzyme_tests = pd.read_csv(tables / "enzyme_family_tests.csv")
    diversity = pd.read_csv(tables / "diversity_tests.csv")

    lines = [
        "# Spatio-temporal profiling of plastic-degrading genetic potential "
        "in Brahmaputra river sediment",
        "",
        f"- Generated (UTC): {utc_now()}",
        "- BioProject: PRJNA944918 (shotgun sediment metagenomes; the SRA `AMPLICON` "
        "label is a known mislabel - see Sharma et al. 2024, DOI 10.3389/fmicb.2024.1426463)",
        "",
        "## Research question",
        "",
        "How does the genetic potential for plastic degradation vary in space (upstream to "
        "downstream, urban versus non-urban) and across seasons in Brahmaputra river sediment, "
        "and which enzyme families and genera drive that variation?",
        "",
    ]
    lines += design_section(summary)
    lines += overview_section(summary, enzyme_raw)
    lines += permanova_section(permanova)
    lines += gradient_section(gradient)
    lines += feature_section(enzyme_tests, diversity)

    fig_dir = Path(args.fig_dir)
    lines += ["## Figures", ""]
    for fig in sorted(fig_dir.glob("*.png")):
        lines.append(f"- `{fig.as_posix()}`")
    lines.append("")
    lines += limitations_section()
    lines += [
        "## Reproducibility",
        "",
        "- Full comparative pipeline: `.venv/bin/python run_comparative.py --threads 12`",
        "- Design table: `refs/sample_metadata.tsv`",
        "- Tables: `outputs/comparative/tables/`, figures: `outputs/comparative/figures/`",
        "",
    ]

    ensure_dir(Path(args.out_md).parent)
    Path(args.out_md).write_text("\n".join(lines), encoding="utf-8")
    print(f"Report written to {args.out_md}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[comparative_report] Error: {exc}", file=sys.stderr)
        sys.exit(1)
