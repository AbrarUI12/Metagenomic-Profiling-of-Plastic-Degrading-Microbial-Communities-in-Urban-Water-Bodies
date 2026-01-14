import argparse
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

from scripts.utils import ensure_dir


def plot_top_genera_bar(raw_matrix, out_path, top_n=20):
    totals = raw_matrix.sum(axis=1).sort_values(ascending=False).head(top_n)
    ensure_dir(Path(out_path).parent)
    plt.figure(figsize=(10, 6))
    sns.barplot(x=totals.values, y=totals.index, color="#3b6ea8")
    plt.xlabel("Total putative plastic-enzyme hits")
    plt.ylabel("Genus")
    plt.title("Top genera by total putative plastic-enzyme hits")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def plot_heatmap(cpm_matrix, raw_matrix, out_path, top_n=30):
    totals = raw_matrix.sum(axis=1).sort_values(ascending=False).head(top_n)
    subset = cpm_matrix.loc[totals.index]
    data = np.log1p(subset)
    ensure_dir(Path(out_path).parent)
    plt.figure(figsize=(12, 8))
    sns.heatmap(data, cmap="viridis")
    plt.title("Genus vs enzyme-family hits (log1p CPM)")
    plt.xlabel("Enzyme family")
    plt.ylabel("Genus")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def choose_k(matrix_scaled, k_min=2, k_max=8):
    inertias = []
    silhouettes = []
    ks = list(range(k_min, k_max + 1))
    for k in ks:
        km = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = km.fit_predict(matrix_scaled)
        inertias.append(km.inertia_)
        if len(set(labels)) > 1 and matrix_scaled.shape[0] > k:
            silhouettes.append(silhouette_score(matrix_scaled, labels))
        else:
            silhouettes.append(np.nan)
    return ks, inertias, silhouettes


def plot_elbow_silhouette(ks, inertias, silhouettes, out_elbow, out_sil):
    ensure_dir(Path(out_elbow).parent)
    plt.figure(figsize=(6, 4))
    plt.plot(ks, inertias, marker="o")
    plt.xlabel("k")
    plt.ylabel("Inertia")
    plt.title("Elbow plot")
    plt.tight_layout()
    plt.savefig(out_elbow, dpi=300)
    plt.close()

    plt.figure(figsize=(6, 4))
    plt.plot(ks, silhouettes, marker="o")
    plt.xlabel("k")
    plt.ylabel("Silhouette score")
    plt.title("Silhouette plot")
    plt.tight_layout()
    plt.savefig(out_sil, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-matrix", default="outputs/tables/genus_enzyme_matrix_raw.csv")
    parser.add_argument("--cpm-matrix", default="outputs/tables/genus_enzyme_matrix_cpm.csv")
    parser.add_argument("--out-clusters", default="outputs/tables/genus_clusters.csv")
    parser.add_argument("--fig-dir", default="outputs/figures")
    parser.add_argument("--k-min", type=int, default=2)
    parser.add_argument("--k-max", type=int, default=8)
    parser.add_argument("--k", type=int, default=0)
    args = parser.parse_args()

    raw_matrix = pd.read_csv(args.raw_matrix, index_col=0)
    cpm_matrix = pd.read_csv(args.cpm_matrix, index_col=0)

    if raw_matrix.empty or raw_matrix.values.sum() == 0:
        ensure_dir(Path(args.fig_dir))
        for name in [
            "top_genera_barplot.png",
            "genus_enzyme_heatmap.png",
            "pca_kmeans.png",
            "kmeans_elbow.png",
            "kmeans_silhouette.png",
        ]:
            plt.figure(figsize=(6, 4))
            plt.text(0.5, 0.5, "No hits to plot", ha="center", va="center")
            plt.axis("off")
            plt.tight_layout()
            plt.savefig(str(Path(args.fig_dir) / name), dpi=300)
            plt.close()
        ensure_dir(Path(args.out_clusters).parent)
        pd.DataFrame(columns=["genus", "cluster_id", "total_hits", "top_families"]).to_csv(
            args.out_clusters, index=False
        )
        return

    plot_top_genera_bar(raw_matrix, str(Path(args.fig_dir) / "top_genera_barplot.png"))
    plot_heatmap(
        cpm_matrix,
        raw_matrix,
        str(Path(args.fig_dir) / "genus_enzyme_heatmap.png"),
    )

    data = np.log1p(cpm_matrix)
    scaler = StandardScaler()
    matrix_scaled = scaler.fit_transform(data.values)

    ks, inertias, silhouettes = choose_k(matrix_scaled, args.k_min, args.k_max)
    plot_elbow_silhouette(
        ks,
        inertias,
        silhouettes,
        str(Path(args.fig_dir) / "kmeans_elbow.png"),
        str(Path(args.fig_dir) / "kmeans_silhouette.png"),
    )

    if args.k and args.k > 1:
        k_final = args.k
    else:
        if all(np.isnan(silhouettes)):
            k_final = 3
        else:
            k_final = ks[int(np.nanargmax(silhouettes))]

    km = KMeans(n_clusters=k_final, random_state=42, n_init=10)
    labels = km.fit_predict(matrix_scaled)

    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(matrix_scaled)
    pca_df = pd.DataFrame(
        {"PC1": coords[:, 0], "PC2": coords[:, 1], "cluster": labels, "genus": raw_matrix.index}
    )

    plt.figure(figsize=(7, 6))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="cluster", palette="tab10")
    plt.title("PCA of genus enzyme profiles (colored by k-means)")
    plt.tight_layout()
    plt.savefig(str(Path(args.fig_dir) / "pca_kmeans.png"), dpi=300)
    plt.close()

    totals = raw_matrix.sum(axis=1)
    top_families = (
        raw_matrix.apply(lambda row: row.sort_values(ascending=False).head(3).index.tolist(), axis=1)
    )

    out_df = pd.DataFrame(
        {
            "genus": raw_matrix.index,
            "cluster_id": labels,
            "total_hits": totals.values,
            "top_families": [";".join(fams) for fams in top_families.values],
        }
    )
    ensure_dir(Path(args.out_clusters).parent)
    out_df.to_csv(args.out_clusters, index=False)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[cluster] Error: {exc}", file=sys.stderr)
        sys.exit(1)
