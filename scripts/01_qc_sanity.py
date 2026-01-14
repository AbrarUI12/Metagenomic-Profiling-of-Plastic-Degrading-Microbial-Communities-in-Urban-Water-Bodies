import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts.utils import ensure_dir, write_json


def fastq_stats(path, sample_reads=10000):
    total_reads = 0
    total_bases = 0
    min_len = None
    max_len = 0
    length_counts = {}
    sample_checked = 0
    sample_mismatch = 0

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline().strip()

            if not qual:
                break

            total_reads += 1
            seq_len = len(seq)
            total_bases += seq_len
            length_counts[seq_len] = length_counts.get(seq_len, 0) + 1
            if min_len is None or seq_len < min_len:
                min_len = seq_len
            if seq_len > max_len:
                max_len = seq_len

            if sample_checked < sample_reads:
                sample_checked += 1
                if len(qual) != seq_len:
                    sample_mismatch += 1

    mean_len = (total_bases / total_reads) if total_reads else 0
    return {
        "path": str(path),
        "total_reads": total_reads,
        "total_bases": total_bases,
        "min_len": min_len,
        "max_len": max_len,
        "mean_len": mean_len,
        "length_counts": length_counts,
        "qc_sample_reads": sample_checked,
        "qc_length_mismatch_reads": sample_mismatch,
    }


def write_markdown(out_md, stats):
    lines = []
    lines.append("# FASTQ sanity summary\n")
    lines.append(f"- Generated (UTC): {datetime.utcnow().isoformat()}Z\n")
    for label, st in stats.items():
        lines.append(f"\n## {label}\n")
        lines.append(f"- Path: {st['path']}\n")
        lines.append(f"- Total reads: {st['total_reads']}\n")
        lines.append(f"- Total bases: {st['total_bases']}\n")
        lines.append(f"- Read length min/max/mean: {st['min_len']} / {st['max_len']} / {st['mean_len']:.2f}\n")
        lines.append(
            f"- QC sample (reads): {st['qc_sample_reads']}, "
            f"length mismatches: {st['qc_length_mismatch_reads']}\n"
        )
        top_lengths = sorted(
            st["length_counts"].items(), key=lambda x: x[1], reverse=True
        )[:5]
        lines.append("- Top read lengths (len:count): " + ", ".join(f"{l}:{c}" for l, c in top_lengths) + "\n")

    ensure_dir(Path(out_md).parent)
    with open(out_md, "w", encoding="utf-8") as f:
        f.writelines(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq1", required=True)
    parser.add_argument("--fastq2", required=False)
    parser.add_argument("--out-json", default="outputs/qc_summary.json")
    parser.add_argument("--out-md", default="outputs/qc_summary.md")
    parser.add_argument("--qc-sample-reads", type=int, default=10000)
    args = parser.parse_args()

    stats = {}
    stats["read1"] = fastq_stats(args.fastq1, sample_reads=args.qc_sample_reads)
    if args.fastq2:
        stats["read2"] = fastq_stats(args.fastq2, sample_reads=args.qc_sample_reads)

    write_json(args.out_json, stats)
    write_markdown(args.out_md, stats)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[qc] Error: {exc}", file=sys.stderr)
        sys.exit(1)
