import argparse
import json
import platform
import subprocess
import sys
from pathlib import Path

from scripts.utils import ensure_dir, run_cmd, which, write_json


def log_versions(out_path):
    lines = []
    lines.append(f"python_version\t{sys.version.split()[0]}\n")
    lines.append(f"os\t{platform.platform()}\n")
    diamond_path = which("diamond")
    local_diamond = Path("tools/diamond/diamond.exe")
    if not diamond_path and local_diamond.exists():
        diamond_path = str(local_diamond)
    if diamond_path:
        lines.append(f"diamond_path\t{diamond_path}\n")
        try:
            out = subprocess.check_output(
                f"\"{diamond_path}\" version", shell=True, text=True, stderr=subprocess.STDOUT
            ).strip()
            lines.append(f"diamond_version\t{out}\n")
        except Exception:
            pass
    for tool in ["mmseqs", "blastx", "makeblastdb", "conda", "pandoc"]:
        tool_path = which(tool)
        if tool_path:
            lines.append(f"{tool}_path\t{tool_path}\n")
    # Python package versions
    pkg_versions = {}
    for pkg in ["pandas", "numpy", "matplotlib", "seaborn", "sklearn", "networkx"]:
        try:
            mod = __import__(pkg)
            ver = getattr(mod, "__version__", "unknown")
            pkg_versions[pkg] = ver
        except Exception:
            pkg_versions[pkg] = "not_installed"
    for k, v in pkg_versions.items():
        lines.append(f"{k}\t{v}\n")

    ensure_dir(Path(out_path).parent)
    Path(out_path).write_text("".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq1", default="SRR23872596_1.fastq")
    parser.add_argument("--fastq2", default="SRR23872596_2.fastq")
    parser.add_argument("--srr", default="SRR23872596")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--pident", type=float, default=30.0)
    parser.add_argument("--min-aln-len", type=int, default=50)
    parser.add_argument("--max-reads", type=int, default=0)
    parser.add_argument("--k-min", type=int, default=2)
    parser.add_argument("--k-max", type=int, default=8)
    parser.add_argument("--k", type=int, default=0)
    parser.add_argument("--render-pdf", action="store_true")
    parser.add_argument("--skip-conda-install", action="store_true")
    parser.add_argument("--network-metric", default="cosine")
    parser.add_argument("--network-threshold", type=float, default=0.6)
    parser.add_argument("--network-presence-threshold", type=float, default=0.0)
    parser.add_argument("--network-algorithm", default="louvain")
    parser.add_argument("--network-seed", type=int, default=42)
    parser.add_argument("--network-log1p", dest="network_log1p", action="store_true")
    parser.add_argument("--no-network-log1p", dest="network_log1p", action="store_false")
    parser.set_defaults(network_log1p=True)
    args = parser.parse_args()

    # Ensure folder structure
    for d in [
        "data",
        "refs",
        "outputs",
        "notebooks",
        "scripts",
        "logs",
        "reports",
        "outputs/figures",
        "outputs/tables",
        "outputs/logs",
        "outputs/enzyme_hits",
        "outputs/metadata",
        "outputs/tmp",
    ]:
        ensure_dir(d)

    commands_log = "outputs/logs/commands.log"
    versions_log = "outputs/logs/versions.txt"

    log_versions(versions_log)

    run_params = {
        "srr": args.srr,
        "fastq1": args.fastq1,
        "fastq2": args.fastq2,
        "threads": args.threads,
        "evalue": args.evalue,
        "pident": args.pident,
        "min_aln_len": args.min_aln_len,
        "max_reads": args.max_reads,
        "k_min": args.k_min,
        "k_max": args.k_max,
        "k": args.k,
        "skip_conda_install": args.skip_conda_install,
        "network_metric": args.network_metric,
        "network_threshold": args.network_threshold,
        "network_presence_threshold": args.network_presence_threshold,
        "network_algorithm": args.network_algorithm,
        "network_seed": args.network_seed,
        "network_log1p": args.network_log1p,
    }
    write_json("outputs/run_params.json", run_params)

    python = sys.executable

    # Step 00: metadata
    sra_table = "SraRunTable.csv"
    meta_cmd = (
        f"\"{python}\" scripts/00_metadata.py --srr {args.srr} "
        f"--out-md outputs/metadata_summary.md --out-json outputs/metadata/metadata.json "
        f"--out-runinfo outputs/metadata/sra_runinfo.csv"
    )
    if Path(sra_table).exists():
        meta_cmd += f" --sra-run-table {sra_table}"
    run_cmd(meta_cmd, commands_log=commands_log)

    # Step 01: QC sanity
    qc_cmd = (
        f"\"{python}\" scripts/01_qc_sanity.py --fastq1 {args.fastq1} "
        f"--fastq2 {args.fastq2} --out-json outputs/qc_summary.json --out-md outputs/qc_summary.md"
    )
    run_cmd(qc_cmd, commands_log=commands_log)

    # Step 02: reference DB
    ref_cmd = (
        f"\"{python}\" scripts/02_fetch_refs.py --out-fasta refs/plastic_enzymes.fasta "
        f"--out-map refs/enzyme_family_map.tsv --out-taxonomy refs/protein_to_taxonomy.tsv "
        f"--out-readme refs/README_refs.md --raw-dir refs/raw"
    )
    run_cmd(ref_cmd, commands_log=commands_log)

    # Step 03: enzyme detection
    search_cmd = (
        f"\"{python}\" scripts/03_diamond_search.py --fastq1 {args.fastq1} --fastq2 {args.fastq2} "
        f"--fasta-db refs/plastic_enzymes.fasta --family-map refs/enzyme_family_map.tsv "
        f"--out outputs/enzyme_hits/enzyme_hits.tsv --threads {args.threads} "
        f"--evalue {args.evalue} --pident {args.pident} --min-aln-len {args.min_aln_len} "
        f"--max-reads {args.max_reads} --commands-log {commands_log} --method-log outputs/logs/search_method.txt"
    )
    if args.skip_conda_install:
        search_cmd += " --skip-conda-install"
    run_cmd(search_cmd, commands_log=commands_log)

    # Step 04: taxonomy mapping
    tax_cmd = (
        f"\"{python}\" scripts/04_taxonomy_map.py --hits outputs/enzyme_hits/enzyme_hits.tsv "
        f"--taxonomy refs/protein_to_taxonomy.tsv"
    )
    run_cmd(tax_cmd, commands_log=commands_log)

    # Step 05: aggregation
    agg_cmd = (
        f"\"{python}\" scripts/05_aggregate.py --hits outputs/enzyme_hits/enzyme_hits.tsv "
        f"--taxonomy refs/protein_to_taxonomy.tsv --qc outputs/qc_summary.json "
        f"--out-hits outputs/tables/enzyme_hits_with_genus.tsv "
        f"--out-raw outputs/tables/genus_enzyme_matrix_raw.csv "
        f"--out-cpm outputs/tables/genus_enzyme_matrix_cpm.csv"
    )
    run_cmd(agg_cmd, commands_log=commands_log)

    # Step 06: clustering and plots
    cluster_cmd = (
        f"\"{python}\" scripts/06_cluster_and_plots.py "
        f"--raw-matrix outputs/tables/genus_enzyme_matrix_raw.csv "
        f"--cpm-matrix outputs/tables/genus_enzyme_matrix_cpm.csv "
        f"--out-clusters outputs/tables/genus_clusters.csv "
        f"--fig-dir outputs/figures --k-min {args.k_min} --k-max {args.k_max} --k {args.k}"
    )
    run_cmd(cluster_cmd, commands_log=commands_log)

    # Step 07: network modularity analysis
    network_cmd = (
        f"\"{python}\" scripts/07_network_modularity.py "
        f"--raw-matrix outputs/tables/genus_enzyme_matrix_raw.csv "
        f"--cpm-matrix outputs/tables/genus_enzyme_matrix_cpm.csv "
        f"--kmeans outputs/tables/genus_clusters_kmeans.csv "
        f"--metric {args.network_metric} --threshold {args.network_threshold} "
        f"--presence-threshold {args.network_presence_threshold} "
        f"--algorithm {args.network_algorithm} --seed {args.network_seed} "
        f"--commands-log {commands_log}"
    )
    if not args.network_log1p:
        network_cmd += " --no-log1p"
    run_cmd(network_cmd, commands_log=commands_log)

    # Step 99: report
    report_cmd = (
        f"\"{python}\" scripts/99_make_report.py --out-md outputs/report.md "
        f"--metadata-json outputs/metadata/metadata.json --qc-json outputs/qc_summary.json "
        f"--raw-matrix outputs/tables/genus_enzyme_matrix_raw.csv "
        f"--cpm-matrix outputs/tables/genus_enzyme_matrix_cpm.csv "
        f"--clusters outputs/tables/genus_clusters.csv --params outputs/run_params.json "
        f"--method-log outputs/logs/search_method.txt"
    )
    if args.render_pdf:
        report_cmd += " --render-pdf"
    run_cmd(report_cmd, commands_log=commands_log)


if __name__ == "__main__":
    main()
