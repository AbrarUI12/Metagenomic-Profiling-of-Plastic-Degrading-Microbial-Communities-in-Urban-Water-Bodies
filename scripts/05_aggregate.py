import argparse
import json
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

import pandas as pd

from scripts.utils import ensure_dir


def total_reads_from_qc(qc_path):
    if not Path(qc_path).exists():
        return 0
    with open(qc_path, "r", encoding="utf-8") as f:
        qc = json.load(f)
    total = 0
    for key in ["read1", "read2"]:
        if key in qc:
            total += int(qc[key].get("total_reads", 0))
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits", default="outputs/enzyme_hits/enzyme_hits.tsv")
    parser.add_argument("--taxonomy", default="refs/protein_to_taxonomy.tsv")
    parser.add_argument("--qc", default="outputs/qc_summary.json")
    parser.add_argument("--out-hits", default="outputs/tables/enzyme_hits_with_genus.tsv")
    parser.add_argument("--out-raw", default="outputs/tables/genus_enzyme_matrix_raw.csv")
    parser.add_argument("--out-cpm", default="outputs/tables/genus_enzyme_matrix_cpm.csv")
    args = parser.parse_args()

    hits = pd.read_csv(args.hits, sep="\t")
    tax = pd.read_csv(args.taxonomy, sep="\t")
    merged = hits.merge(tax, on="protein_id", how="left")
    merged["genus"] = merged["genus"].fillna("Unassigned")
    ensure_dir(Path(args.out_hits).parent)
    merged.to_csv(args.out_hits, sep="\t", index=False)

    matrix = (
        merged.groupby(["genus", "family"])
        .size()
        .reset_index(name="count")
        .pivot(index="genus", columns="family", values="count")
        .fillna(0)
        .astype(int)
    )
    matrix.to_csv(args.out_raw)

    total_reads = total_reads_from_qc(args.qc)
    if total_reads > 0:
        cpm = matrix / total_reads * 1e6
    else:
        cpm = matrix.copy()
    cpm.to_csv(args.out_cpm)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[aggregate] Error: {exc}", file=sys.stderr)
        sys.exit(1)
