"""Figures for the spatio-temporal comparison.

Every figure is driven by the CSVs written by ``12_comparative_stats.py`` and
``11_build_matrices.py`` - nothing is recomputed here, so the figures and the
reported statistics can never disagree.

Sites are ordered and coloured upstream -> downstream so the spatial axis reads the
same way in every panel.
"""

import argparse
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts.utils import ensure_dir

SEASON_ORDER = ["post_monsoon", "winter", "pre_monsoon"]
SEASON_MARKERS = {"post_monsoon": "o", "winter": "s", "pre_monsoon": "^"}
SEASON_LABELS = {
    "post_monsoon": "Post-monsoon (Sep 2021)",
    "winter": "Winter (Dec 2021-Jan 2022)",
    "pre_monsoon": "Pre-monsoon (Apr 2022)",
}


def site_order(summary):
    """Sites sorted upstream (high latitude) -> downstream (low latitude)."""
    return (
        summary.groupby("site")["lat"].first().sort_values(ascending=False).index.tolist()
    )


def site_colors(sites):
    cmap = plt.get_cmap("viridis")
    return {site: cmap(i / max(len(sites) - 1, 1)) for i, site in enumerate(sites)}


def seasons_present(summary):
    return [s for s in SEASON_ORDER if s in set(summary["season"])]


def plot_pcoa(coords, explained, colors, summary, out_path, title):
    fig, ax = plt.subplots(figsize=(8, 6.5))
    for season in seasons_present(coords):
        sub = coords[coords["season"] == season]
        for site, rows in sub.groupby("site"):
            ax.scatter(
                rows["PCo1"], rows["PCo2"],
                color=colors[site], marker=SEASON_MARKERS.get(season, "o"),
                s=140, edgecolor="black", linewidth=0.6, zorder=3,
            )
    for run, row in coords.iterrows():
        ax.annotate(
            row["site"][:4], (row["PCo1"], row["PCo2"]),
            fontsize=7, xytext=(6, 4), textcoords="offset points", alpha=0.75,
        )
    ax.axhline(0, color="grey", lw=0.5, zorder=1)
    ax.axvline(0, color="grey", lw=0.5, zorder=1)
    ax.set_xlabel(f"PCo1 ({explained[0]:.1%} of variation)")
    ax.set_ylabel(f"PCo2 ({explained[1]:.1%} of variation)")
    ax.set_title(title)

    site_handles = [
        plt.Line2D([], [], marker="o", linestyle="", color=colors[s],
                   markeredgecolor="black", markersize=9, label=s)
        for s in colors
    ]
    season_handles = [
        plt.Line2D([], [], marker=SEASON_MARKERS[s], linestyle="", color="grey",
                   markeredgecolor="black", markersize=9, label=SEASON_LABELS[s])
        for s in seasons_present(coords)
    ]
    # Wrapped title: a single-line version overflows the figure and gets clipped.
    first = ax.legend(handles=site_handles, title="Site\n(upstream -> downstream)",
                      loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8)
    ax.add_artist(first)
    ax.legend(handles=season_handles, title="Season", loc="lower left",
              bbox_to_anchor=(1.02, 0), fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_gradient(summary, colors, out_path):
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for season in seasons_present(summary):
        sub = summary[summary["season"] == season].sort_values("lat", ascending=False)
        ax.plot(sub["lat"], sub["total_hits_cpm"], linestyle="--", alpha=0.5,
                color="grey", zorder=1)
        ax.scatter(sub["lat"], sub["total_hits_cpm"],
                   marker=SEASON_MARKERS.get(season, "o"),
                   c=[colors[s] for s in sub["site"]], s=140,
                   edgecolor="black", linewidth=0.6, zorder=3,
                   label=SEASON_LABELS[season])
    for _, row in summary.iterrows():
        ax.annotate(row["site"][:4], (row["lat"], row["total_hits_cpm"]),
                    fontsize=7, xytext=(5, 5), textcoords="offset points", alpha=0.75)
    # Latitude decreases downstream; invert so the river reads left-to-right.
    ax.invert_xaxis()
    ax.set_xlabel("Latitude (upstream -> downstream)")
    ax.set_ylabel("Plastic-enzyme hits (CPM)")
    ax.set_title("Plastic-degrading potential along the Brahmaputra")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_seasonal_box(summary, enzyme_cpm, out_path):
    families = list(enzyme_cpm.columns)
    seasons = seasons_present(summary)
    n_panels = len(families) + 1
    ncols = min(3, n_panels)
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 3.8 * nrows), squeeze=False)

    panels = [("Total plastic-enzyme potential", summary["total_hits_cpm"])]
    panels += [(f, enzyme_cpm[f]) for f in families]

    for idx, (title, values) in enumerate(panels):
        ax = axes[idx // ncols][idx % ncols]
        data = [values.loc[summary.index[summary["season"] == s]].to_numpy() for s in seasons]
        ax.boxplot(data, tick_labels=[s.replace("_", "\n") for s in seasons], showfliers=False)
        for i, arr in enumerate(data, start=1):
            jitter = np.random.default_rng(42).normal(0, 0.045, len(arr))
            ax.scatter(np.full(len(arr), i) + jitter, arr, s=32, alpha=0.85,
                       edgecolor="black", linewidth=0.4, zorder=3)
        ax.set_title(title, fontsize=10)
        ax.set_ylabel("CPM")
    for j in range(len(panels), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.suptitle("Seasonal variation in plastic-degrading potential", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_heatmap(enzyme_cpm, summary, out_path):
    order = summary.sort_values(["lat", "season"], ascending=[False, True]).index
    data = enzyme_cpm.loc[order]
    labels = [
        f"{summary.loc[r, 'site'][:9]} | {summary.loc[r, 'season'][:4]}" for r in order
    ]
    fig, ax = plt.subplots(figsize=(1.5 + 1.35 * len(data.columns), 0.42 * len(order) + 2.2))
    im = ax.imshow(data.to_numpy(), aspect="auto", cmap="magma")
    ax.set_xticks(range(len(data.columns)))
    ax.set_xticklabels(data.columns, rotation=40, ha="right")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(labels, fontsize=8)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            value = data.iat[i, j]
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=6.5,
                    color="white" if value < data.to_numpy().max() * 0.6 else "black")
    fig.colorbar(im, ax=ax, label="CPM")
    ax.set_title("Enzyme-family potential per sample\n(upstream -> downstream)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_composition(enzyme_cpm, summary, out_path):
    order = summary.sort_values(["lat", "season"], ascending=[False, True]).index
    data = enzyme_cpm.loc[order]
    totals = data.sum(axis=1).replace(0, np.nan)
    fractions = data.div(totals, axis=0).fillna(0)
    labels = [
        f"{summary.loc[r, 'site'][:9]} | {summary.loc[r, 'season'][:4]}" for r in order
    ]
    fig, ax = plt.subplots(figsize=(11, 0.45 * len(order) + 2.5))
    left = np.zeros(len(order))
    cmap = plt.get_cmap("tab10")
    for i, family in enumerate(data.columns):
        ax.barh(range(len(order)), fractions[family], left=left,
                color=cmap(i % 10), edgecolor="white", label=family)
        left += fractions[family].to_numpy()
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Fraction of plastic-enzyme hits")
    ax.set_title("Enzyme-family composition per sample (upstream -> downstream)")
    ax.legend(fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables-dir", default="outputs/comparative/tables")
    parser.add_argument("--fig-dir", default="outputs/comparative/figures")
    args = parser.parse_args()

    tables = Path(args.tables_dir)
    fig_dir = Path(args.fig_dir)
    ensure_dir(fig_dir)

    summary = pd.read_csv(tables / "sample_summary.tsv", sep="\t", index_col="run_accession")
    enzyme_cpm = pd.read_csv(tables / "sample_enzyme_matrix_cpm.csv", index_col="run_accession")

    sites = site_order(summary)
    colors = site_colors(sites)

    for name, title in [
        ("enzyme", "Enzyme-family composition (Bray-Curtis PCoA)"),
        ("genus", "Genus composition of enzyme hits (Bray-Curtis PCoA)"),
    ]:
        coords_path = tables / f"pcoa_{name}_coords.csv"
        if not coords_path.exists():
            continue
        coords = pd.read_csv(coords_path, index_col="run_accession")
        explained = pd.read_csv(tables / f"pcoa_{name}_explained.csv")["proportion_explained"]
        # A 2-D ordination needs at least two axes, i.e. three or more samples.
        if "PCo2" not in coords.columns or len(explained) < 2:
            print(f"  skipping pcoa_{name}.png (too few samples for a 2-D ordination)")
            continue
        plot_pcoa(coords, explained.tolist(), colors, summary,
                  fig_dir / f"pcoa_{name}.png", title)
        print(f"  wrote pcoa_{name}.png")

    plot_gradient(summary, colors, fig_dir / "gradient_latitude.png")
    plot_seasonal_box(summary, enzyme_cpm, fig_dir / "seasonal_boxplots.png")
    plot_heatmap(enzyme_cpm, summary, fig_dir / "enzyme_heatmap.png")
    plot_composition(enzyme_cpm, summary, fig_dir / "enzyme_composition.png")
    print(f"Figures written to {fig_dir}/")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[comparative_plots] Error: {exc}", file=sys.stderr)
        sys.exit(1)
