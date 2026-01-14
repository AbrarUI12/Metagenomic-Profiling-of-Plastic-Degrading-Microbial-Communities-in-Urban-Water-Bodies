import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

import pandas as pd

from scripts.utils import ensure_dir, log_command, which


TOPIC_DESCRIPTION = (
    "Using shotgun metagenomic reads from a polluted freshwater environment, screen "
    "for homologs of known plastic-degrading enzymes (e.g., PETase/cutinase-like "
    "hydrolases, MHETase-like, alkane monooxygenase AlkB, relevant esterases/lipases "
    "where appropriate). Assign each enzyme hit to a likely microbial genus by "
    "best-hit taxonomy (homology), quantify enzyme-hit abundance per genus, then "
    "perform k-means clustering of genera based on enzyme-hit profiles to identify "
    "'high / moderate / low plastic-degradation potential' groups. Report results "
    "with conservative, thesis-safe language emphasizing genetic potential and "
    "limitations."
)


def load_json(path):
    if not Path(path).exists():
        return {}
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def fmt_number(n):
    try:
        return f"{int(n):,}"
    except Exception:
        return str(n)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-md", default="outputs/report.md")
    parser.add_argument("--metadata-json", default="outputs/metadata/metadata.json")
    parser.add_argument("--qc-json", default="outputs/qc_summary.json")
    parser.add_argument("--raw-matrix", default="outputs/tables/genus_enzyme_matrix_raw.csv")
    parser.add_argument("--cpm-matrix", default="outputs/tables/genus_enzyme_matrix_cpm.csv")
    parser.add_argument("--clusters", default="outputs/tables/genus_clusters.csv")
    parser.add_argument("--params", default="outputs/run_params.json")
    parser.add_argument("--method-log", default="outputs/logs/search_method.txt")
    parser.add_argument("--render-pdf", action="store_true")
    parser.add_argument("--commands-log", default="outputs/logs/commands.log")
    args = parser.parse_args()

    metadata = load_json(args.metadata_json)
    qc = load_json(args.qc_json)
    params = load_json(args.params)
    method_log = Path(args.method_log).read_text(encoding="utf-8").strip() if Path(args.method_log).exists() else ""

    raw = pd.read_csv(args.raw_matrix, index_col=0)
    cpm = pd.read_csv(args.cpm_matrix, index_col=0)
    clusters = pd.read_csv(args.clusters)

    total_hits = int(raw.values.sum())
    genera_with_hits = int((raw.sum(axis=1) > 0).sum())

    top_genera = raw.sum(axis=1).sort_values(ascending=False).head(10)
    top_families = raw.sum(axis=0).sort_values(ascending=False).head(10)

    read1 = qc.get("read1", {})
    read2 = qc.get("read2", {})
    total_reads = int(read1.get("total_reads", 0)) + int(read2.get("total_reads", 0))

    lines = []
    lines.append("# Plastic-degrading genetic potential profiling report\n")
    lines.append(f"- Generated (UTC): {datetime.utcnow().isoformat()}Z\n")
    lines.append("\n## Research topic\n")
    lines.append(f"{TOPIC_DESCRIPTION}\n")

    lines.append("\n## Data description\n")
    lines.append(f"- Accession: {metadata.get('srr', 'SRR23872596')}\n")
    if total_reads:
        lines.append(f"- Total reads (R1+R2): {fmt_number(total_reads)}\n")
    if read1:
        lines.append(
            f"- Read1 length (min/mean/max): {read1.get('min_len')} / "
            f"{read1.get('mean_len', 0):.2f} / {read1.get('max_len')}\n"
        )
    if read2:
        lines.append(
            f"- Read2 length (min/mean/max): {read2.get('min_len')} / "
            f"{read2.get('mean_len', 0):.2f} / {read2.get('max_len')}\n"
        )
    summary = metadata.get("summary", {})
    if summary:
        lines.append("- Key metadata fields:\n")
        for k, v in summary.items():
            lines.append(f"  - {k}: {v}\n")

    lines.append("\n## Methods (homology-based inference)\n")
    lines.append(
        "- Reference database: UniProtKB stream queries for PETase-like, MHETase-like, "
        "cutinase-like, AlkB, and polyesterase-like enzymes (see `refs/README_refs.md`).\n"
    )
    lines.append(
        "- Enzyme detection: translated search of reads vs protein database; top hit per read retained.\n"
    )
    if params:
        lines.append(
            "- Thresholds: e-value <= {evalue}, percent identity >= {pident}, min alignment length >= {min_aln_len}.\n".format(
                evalue=params.get("evalue"),
                pident=params.get("pident"),
                min_aln_len=params.get("min_aln_len"),
            )
        )
        if params.get("max_reads", 0):
            lines.append(
                f"- Subsampling: analysis limited to first {params.get('max_reads')} reads per FASTQ for this run.\n"
            )
    if method_log:
        lines.append(f"- Search method log: `{args.method_log}`\n")
    lines.append(
        "- Taxonomy assignment: best-hit mapping from UniProt protein accessions to genus (homology-based).\n"
    )
    lines.append(
        "- Clustering: k-means on log1p CPM enzyme-family profiles; k selected by silhouette.\n"
    )

    lines.append("\n## Results (putative genetic potential)\n")
    lines.append(
        f"- Total putative enzyme hits: {fmt_number(total_hits)} across {fmt_number(genera_with_hits)} genera.\n"
    )
    lines.append("- Top genera by total putative plastic-enzyme hits:\n")
    for genus, count in top_genera.items():
        lines.append(f"  - {genus}: {fmt_number(count)}\n")
    lines.append("- Top enzyme families by hit count:\n")
    for fam, count in top_families.items():
        lines.append(f"  - {fam}: {fmt_number(count)}\n")
    lines.append("\nCluster summary (putative plastic-degradation potential groups):\n")
    cluster_counts = clusters["cluster_id"].value_counts().sort_index()
    for cid, count in cluster_counts.items():
        lines.append(f"- Cluster {cid}: {fmt_number(count)} genera\n")

    lines.append("\n## Figures\n")
    lines.append("- `outputs/figures/top_genera_barplot.png`\n")
    lines.append("- `outputs/figures/genus_enzyme_heatmap.png`\n")
    lines.append("- `outputs/figures/pca_kmeans.png`\n")
    lines.append("- `outputs/figures/kmeans_elbow.png`\n")
    lines.append("- `outputs/figures/kmeans_silhouette.png`\n")

    lines.append("\n## Limitations (thesis-safe)\n")
    lines.append(
        "- These are putative, candidate plastic-degrading genes identified by homology-based inference; "
        "metagenomic presence does not imply expression or active degradation.\n"
    )
    lines.append(
        "- Best-hit taxonomy can misassign genes across taxa, especially for conserved hydrolases.\n"
    )
    lines.append(
        "- The database is curated from public annotations and may include broad enzyme families; "
        "functional specificity requires experimental validation.\n"
    )
    lines.append(
        "- Results describe genetic potential in bulk sediment metagenome and do not prove plastisphere association.\n"
    )

    lines.append("\n## Reproducibility\n")
    lines.append("- Run end-to-end pipeline: `python run_pipeline.py`\n")
    lines.append("- Outputs are saved under `outputs/` with logs in `outputs/logs/`.\n")

    ensure_dir(Path(args.out_md).parent)
    Path(args.out_md).write_text("".join(lines), encoding="utf-8")

    if args.render_pdf and which("pandoc"):
        pdf_path = Path(args.out_md).with_suffix(".pdf")
        cmd = f"pandoc \"{args.out_md}\" -o \"{pdf_path}\""
        log_command(args.commands_log, cmd)
        import subprocess

        subprocess.run(cmd, shell=True, check=False)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[report] Error: {exc}", file=sys.stderr)
        sys.exit(1)
