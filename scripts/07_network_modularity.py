import argparse
import json
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns

from scripts.utils import ensure_dir, log_command, write_json


def load_matrix(path):
    return pd.read_csv(path, index_col=0)


def log1p_if_needed(data, use_log1p):
    if use_log1p:
        return np.log1p(data)
    return data


def cosine_similarity_matrix(data):
    X = data.values.astype(float)
    norms = np.linalg.norm(X, axis=1)
    norms[norms == 0] = 1.0
    Xn = X / norms[:, None]
    sim = Xn @ Xn.T
    np.fill_diagonal(sim, 1.0)
    return pd.DataFrame(sim, index=data.index, columns=data.index)


def corr_similarity_matrix(data, method):
    sim = data.T.corr(method=method)
    sim = sim.fillna(0.0)
    np.fill_diagonal(sim.values, 1.0)
    return sim


def jaccard_similarity_matrix(data, presence_threshold):
    binary = (data.values > presence_threshold).astype(int)
    intersection = binary @ binary.T
    row_sums = binary.sum(axis=1)
    union = row_sums[:, None] + row_sums[None, :] - intersection
    with np.errstate(divide="ignore", invalid="ignore"):
        sim = np.where(union > 0, intersection / union, 0.0)
    np.fill_diagonal(sim, 1.0)
    return pd.DataFrame(sim, index=data.index, columns=data.index)


def build_graph(sim_df, threshold):
    G = nx.Graph()
    genera = list(sim_df.index)
    G.add_nodes_from(genera)
    for i, g1 in enumerate(genera):
        for j in range(i + 1, len(genera)):
            g2 = genera[j]
            w = float(sim_df.iat[i, j])
            if w > threshold:
                G.add_edge(g1, g2, weight=w)
    return G


def detect_communities(G, algorithm, seed):
    if G.number_of_edges() == 0:
        return [{n} for n in G.nodes()], "none"

    if algorithm == "louvain":
        try:
            from networkx.algorithms.community import louvain_communities

            communities = louvain_communities(G, weight="weight", seed=seed)
            return communities, "louvain"
        except Exception:
            pass

    from networkx.algorithms.community import greedy_modularity_communities

    communities = greedy_modularity_communities(G, weight="weight")
    return communities, "greedy_modularity"


def assign_module_ids(communities):
    def sort_key(comm):
        return (-len(comm), sorted(comm)[0])

    ordered = sorted(communities, key=sort_key)
    module_map = {}
    for idx, comm in enumerate(ordered):
        for genus in comm:
            module_map[genus] = idx
    return module_map, ordered


def top_families_from_row(row, top_n=3):
    row = row.sort_values(ascending=False)
    row = row[row > 0]
    if row.empty:
        return ""
    return ";".join(row.head(top_n).index.tolist())


def compute_module_summary(raw_matrix, module_map):
    df = raw_matrix.copy()
    df["module_id"] = df.index.map(module_map)
    summaries = []
    for module_id, group in df.groupby("module_id"):
        counts = group.drop(columns=["module_id"]).sum(axis=0)
        top_fams = counts.sort_values(ascending=False)
        top_fams = top_fams[top_fams > 0].head(3).index.tolist()
        summaries.append(
            {
                "module_id": module_id,
                "num_genera": int(group.shape[0]),
                "total_hits": int(group.drop(columns=["module_id"]).values.sum()),
                "top_families": ";".join(top_fams),
            }
        )
    return pd.DataFrame(summaries).sort_values("module_id")


def save_network_tables(G, module_map, raw_matrix, edges_path, nodes_path):
    ensure_dir(Path(edges_path).parent)
    edges = []
    for u, v, data in G.edges(data=True):
        edges.append([u, v, data.get("weight", 0.0)])
    pd.DataFrame(edges, columns=["genus1", "genus2", "similarity"]).to_csv(
        edges_path, sep="\t", index=False
    )

    degree = dict(G.degree(weight="weight"))
    totals = raw_matrix.sum(axis=1)
    nodes = []
    for genus in raw_matrix.index:
        nodes.append(
            [
                genus,
                module_map.get(genus, -1),
                float(degree.get(genus, 0.0)),
                int(totals.get(genus, 0)),
            ]
        )
    pd.DataFrame(
        nodes, columns=["genus", "module_id", "degree", "total_enzyme_hits"]
    ).to_csv(nodes_path, sep="\t", index=False)


def plot_network(G, module_map, raw_matrix, out_path, seed):
    ensure_dir(Path(out_path).parent)
    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(G, seed=seed, weight="weight")

    modules = sorted(set(module_map.values()))
    palette = sns.color_palette("tab20", max(1, len(modules)))
    color_map = {m: palette[i % len(palette)] for i, m in enumerate(modules)}

    totals = raw_matrix.sum(axis=1)
    node_sizes = []
    node_colors = []
    for node in G.nodes():
        node_sizes.append(60 + 40 * np.log1p(totals.get(node, 0)))
        node_colors.append(color_map.get(module_map.get(node, -1), (0.5, 0.5, 0.5)))

    weights = [G[u][v].get("weight", 0.0) for u, v in G.edges()]
    if weights:
        max_w = max(weights)
        widths = [0.5 + 2.0 * (w / max_w) for w in weights]
    else:
        widths = []

    nx.draw_networkx_edges(G, pos, alpha=0.4, width=widths, edge_color="#999999")
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors)

    labels = {}
    if G.number_of_nodes() <= 50:
        labels = {n: n for n in G.nodes()}
    else:
        top_nodes = totals.sort_values(ascending=False).head(20).index
        labels = {n: n for n in top_nodes if n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=7)

    plt.title("Genus similarity network (modules colored)")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def plot_heatmap_by_module(cpm_matrix, module_map, out_path, use_log1p):
    ensure_dir(Path(out_path).parent)
    order = sorted(cpm_matrix.index, key=lambda g: (module_map.get(g, 9999), g))
    data = cpm_matrix.loc[order]
    data = log1p_if_needed(data, use_log1p)
    plt.figure(figsize=(12, 8))
    sns.heatmap(data, cmap="viridis")
    plt.title("Genus vs enzyme families ordered by module")
    plt.xlabel("Enzyme family")
    plt.ylabel("Genus (ordered by module)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def plot_module_bars(summary_df, out_path):
    ensure_dir(Path(out_path).parent)
    summary_df = summary_df.sort_values("module_id")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sns.barplot(
        x="module_id",
        y="num_genera",
        data=summary_df,
        ax=axes[0],
        color="#4c78a8",
    )
    axes[0].set_title("Number of genera per module")
    axes[0].set_xlabel("Module")
    axes[0].set_ylabel("Genera")

    sns.barplot(
        x="module_id",
        y="total_hits",
        data=summary_df,
        ax=axes[1],
        color="#f58518",
    )
    axes[1].set_title("Total enzyme hits per module")
    axes[1].set_xlabel("Module")
    axes[1].set_ylabel("Total hits")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def plot_kmeans_vs_modules(modules_df, kmeans_path, out_path, out_table):
    if not kmeans_path or not Path(kmeans_path).exists():
        return None

    kmeans = pd.read_csv(kmeans_path)
    if "genus" not in kmeans.columns or "cluster_id" not in kmeans.columns:
        return None

    merged = modules_df.merge(kmeans, on="genus", how="left")
    merged = merged.dropna(subset=["cluster_id"])
    if merged.empty:
        return None

    table = (
        merged.groupby(["module_id", "cluster_id"])
        .size()
        .reset_index(name="count")
        .pivot(index="module_id", columns="cluster_id", values="count")
        .fillna(0)
        .astype(int)
    )
    ensure_dir(Path(out_table).parent)
    table.to_csv(out_table)

    ensure_dir(Path(out_path).parent)
    plt.figure(figsize=(8, 6))
    sns.heatmap(table, annot=True, fmt="d", cmap="Blues")
    plt.title("K-means clusters vs network modules")
    plt.xlabel("K-means cluster")
    plt.ylabel("Module")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    return table


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-matrix", default="outputs/tables/genus_enzyme_matrix_raw.csv")
    parser.add_argument("--cpm-matrix", default="outputs/tables/genus_enzyme_matrix_cpm.csv")
    parser.add_argument("--kmeans", default="outputs/tables/genus_clusters_kmeans.csv")
    parser.add_argument("--metric", default="cosine")
    parser.add_argument("--threshold", type=float, default=0.6)
    parser.add_argument("--presence-threshold", type=float, default=0.0)
    parser.add_argument("--algorithm", default="louvain")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--no-log1p", action="store_true")
    parser.add_argument("--network-dir", default="outputs/network")
    parser.add_argument("--fig-dir", default="outputs/figures/network")
    parser.add_argument("--modules-out", default="outputs/tables/genus_modules_modularity.csv")
    parser.add_argument("--module-summary-out", default="outputs/tables/module_summary_modularity.csv")
    parser.add_argument("--edges-out", default="outputs/network/genus_similarity_edges.tsv")
    parser.add_argument("--nodes-out", default="outputs/network/genus_similarity_nodes.tsv")
    parser.add_argument("--comparison-out", default="outputs/tables/kmeans_vs_modules.csv")
    parser.add_argument("--commands-log", default="outputs/logs/commands.log")
    parser.add_argument("--params-out", default="outputs/logs/network_modularity_params.json")
    args = parser.parse_args()

    use_log1p = not args.no_log1p
    np.random.seed(args.seed)

    if args.commands_log:
        log_command(
            args.commands_log,
            f"network_modularity metric={args.metric} threshold={args.threshold} "
            f"presence_threshold={args.presence_threshold} algorithm={args.algorithm} "
            f"seed={args.seed} log1p={use_log1p}",
        )

    raw_matrix = load_matrix(args.raw_matrix)
    cpm_matrix = load_matrix(args.cpm_matrix)

    data = cpm_matrix.copy()
    data = log1p_if_needed(data, use_log1p)

    metric = args.metric.lower()
    if metric == "cosine":
        sim_df = cosine_similarity_matrix(data)
    elif metric in ["pearson", "spearman"]:
        sim_df = corr_similarity_matrix(data, metric)
    elif metric == "jaccard":
        sim_df = jaccard_similarity_matrix(data, args.presence_threshold)
    else:
        raise ValueError("Unsupported metric: " + args.metric)

    G = build_graph(sim_df, args.threshold)
    communities, algo_used = detect_communities(G, args.algorithm, args.seed)
    module_map, ordered = assign_module_ids(communities)

    if G.number_of_edges() == 0:
        modularity_score = 0.0
    else:
        modularity_score = nx.algorithms.community.modularity(G, ordered, weight="weight")

    totals = raw_matrix.sum(axis=1)
    degree = dict(G.degree(weight="weight"))
    modules_df = pd.DataFrame(
        {
            "genus": raw_matrix.index,
            "module_id": [module_map.get(g, -1) for g in raw_matrix.index],
            "degree": [float(degree.get(g, 0.0)) for g in raw_matrix.index],
            "total_enzyme_hits": [int(totals.get(g, 0)) for g in raw_matrix.index],
            "dominant_enzyme_families": [
                top_families_from_row(raw_matrix.loc[g]) for g in raw_matrix.index
            ],
        }
    )
    ensure_dir(Path(args.modules_out).parent)
    modules_df.to_csv(args.modules_out, index=False)

    module_summary = compute_module_summary(raw_matrix, module_map)
    ensure_dir(Path(args.module_summary_out).parent)
    module_summary.to_csv(args.module_summary_out, index=False)

    save_network_tables(G, module_map, raw_matrix, args.edges_out, args.nodes_out)

    plot_network(G, module_map, raw_matrix, str(Path(args.fig_dir) / "genus_network_modules.png"), args.seed)
    plot_heatmap_by_module(
        cpm_matrix,
        module_map,
        str(Path(args.fig_dir) / "genus_enzyme_heatmap_by_module.png"),
        use_log1p,
    )
    plot_module_bars(module_summary, str(Path(args.fig_dir) / "module_barplots.png"))

    kmeans_path = args.kmeans
    if not Path(kmeans_path).exists():
        alt = Path("outputs/tables/genus_clusters.csv")
        if alt.exists():
            kmeans_path = str(alt)
        else:
            kmeans_path = ""
    table = plot_kmeans_vs_modules(
        modules_df,
        kmeans_path,
        str(Path(args.fig_dir) / "kmeans_vs_modules_heatmap.png"),
        args.comparison_out,
    )

    ari = None
    if table is not None and kmeans_path:
        try:
            from sklearn.metrics import adjusted_rand_score

            merged = modules_df.merge(
                pd.read_csv(kmeans_path), on="genus", how="left"
            ).dropna(subset=["cluster_id"])
            if not merged.empty:
                ari = adjusted_rand_score(merged["module_id"], merged["cluster_id"])
        except Exception:
            ari = None

    params = {
        "metric": args.metric,
        "threshold": args.threshold,
        "presence_threshold": args.presence_threshold,
        "algorithm_requested": args.algorithm,
        "algorithm_used": algo_used,
        "seed": args.seed,
        "log1p": use_log1p,
        "num_nodes": int(G.number_of_nodes()),
        "num_edges": int(G.number_of_edges()),
        "num_modules": int(len(set(module_map.values()))),
        "modularity_Q": float(modularity_score),
        "adjusted_rand_index": None if ari is None else float(ari),
    }
    write_json(args.params_out, params)


if __name__ == "__main__":
    main()

