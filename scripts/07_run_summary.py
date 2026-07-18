# 07_run_summary.py
import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

# Allow "from scripts.utils import ensure_dir, write_json" like your other files
sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts.utils import ensure_dir, write_json


def load_json(path: str) -> dict:
    p = Path(path)
    if not p.exists():
        return {}
    with p.open("r", encoding="utf-8") as f:
        return json.load(f)


def read_family_map(path: str) -> dict:
    """protein_id -> family"""
    p = Path(path)
    if not p.exists():
        return {}
    mapping = {}
    with p.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pid = (row.get("protein_id") or "").strip()
            fam = (row.get("family") or "").strip()
            if pid:
                mapping[pid] = fam or "unknown"
    return mapping


def read_hits(hits_tsv: str) -> list[dict]:
    p = Path(hits_tsv)
    if not p.exists():
        return []
    with p.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def qc_totals(qc: dict) -> dict:
    """
    Returns totals for reads/bases and per-readfile totals if available.
    """
    read1 = qc.get("read1", {}) if isinstance(qc, dict) else {}
    read2 = qc.get("read2", {}) if isinstance(qc, dict) else {}

    def safe_int(x, default=0):
        try:
            return int(float(x))
        except Exception:
            return default

    total_reads_r1 = safe_int(read1.get("total_reads", 0))
    total_reads_r2 = safe_int(read2.get("total_reads", 0))
    total_bases_r1 = safe_int(read1.get("total_bases", 0))
    total_bases_r2 = safe_int(read2.get("total_bases", 0))

    return {
        "read1_total_reads": total_reads_r1,
        "read2_total_reads": total_reads_r2,
        "total_reads": total_reads_r1 + total_reads_r2,
        "read1_total_bases": total_bases_r1,
        "read2_total_bases": total_bases_r2,
        "total_bases": total_bases_r1 + total_bases_r2,
        "read1_path": read1.get("path", ""),
        "read2_path": read2.get("path", ""),
    }


def estimate_reads_used(qc_summary: dict, run_params: dict) -> dict:
    """
    If run_params includes max_reads (subsampling per FASTQ in your search script),
    estimate reads used as sum(min(total_reads_file, max_reads)) across read1/read2.

    If max_reads missing or 0 => assume all reads used.
    """
    totals = qc_totals(qc_summary)
    max_reads = 0
    try:
        max_reads = int(run_params.get("max_reads", 0))
    except Exception:
        max_reads = 0

    r1_total = totals["read1_total_reads"]
    r2_total = totals["read2_total_reads"]

    if max_reads and max_reads > 0:
        used_r1 = min(r1_total, max_reads) if r1_total else 0
        used_r2 = min(r2_total, max_reads) if r2_total else 0
        note = f"Subsampled to max_reads={max_reads} per FASTQ (estimated from qc totals)."
    else:
        used_r1 = r1_total
        used_r2 = r2_total
        note = "No subsampling detected (max_reads missing/0). Assuming all reads used."

    return {
        "max_reads_param": max_reads,
        "estimated_used_reads_read1": used_r1,
        "estimated_used_reads_read2": used_r2,
        "estimated_used_reads_total": used_r1 + used_r2,
        "note": note,
    }


def normalize_protein_id(pid: str) -> str:
    """
    Handle cases like sp|QXXXX|NAME where the accession is in the middle.
    Your pipeline already normalizes in 03_diamond_search.py, but this keeps it safe.
    """
    if not pid:
        return ""
    pid = pid.strip()
    if "|" in pid:
        parts = pid.split("|")
        if len(parts) >= 2 and parts[1]:
            return parts[1].strip()
    return pid


def main():
    ap = argparse.ArgumentParser(description="Summarize data usage and plastic-enzyme hit counts.")
    ap.add_argument("--qc", default="outputs/qc_summary.json")
    ap.add_argument("--params", default="outputs/run_params.json")
    ap.add_argument("--hits", default="outputs/enzyme_hits/enzyme_hits.tsv")
    ap.add_argument("--family-map", default="refs/enzyme_family_map.tsv")
    ap.add_argument("--out-json", default="outputs/summary/run_summary.json")
    ap.add_argument("--out-md", default="outputs/summary/run_summary.md")
    args = ap.parse_args()

    qc = load_json(args.qc)
    params = load_json(args.params)
    hits = read_hits(args.hits)
    fam_map = read_family_map(args.family_map)

    totals = qc_totals(qc)
    used = estimate_reads_used(qc, params)

    # Hit-based metrics
    total_hit_rows = len(hits)
    unique_read_ids = set()
    unique_protein_ids = set()
    families_hit = {}

    for row in hits:
        rid = (row.get("read_id") or "").strip()
        pid = normalize_protein_id(row.get("protein_id", ""))
        if rid:
            unique_read_ids.add(rid)
        if pid:
            unique_protein_ids.add(pid)
            fam = row.get("family") or fam_map.get(pid, "unknown")
            fam = (fam or "unknown").strip()
            families_hit[fam] = families_hit.get(fam, 0) + 1

    unique_families = sorted([k for k in families_hit.keys() if k])

    # Some people like “genes” = unique proteins; report both flavors
    summary = {
        "generated_utc": datetime.now(timezone.utc).replace(tzinfo=None).isoformat() + "Z",
        "inputs": {
            "qc_path": args.qc,
            "params_path": args.params,
            "hits_path": args.hits,
            "family_map_path": args.family_map,
        },
        "sequencing_totals_from_qc": totals,
        "estimated_reads_used": used,
        "plastic_enzyme_hits": {
            "total_hit_rows": total_hit_rows,
            "unique_reads_with_hits": len(unique_read_ids),
            "unique_plastic_related_proteins_found": len(unique_protein_ids),
            "unique_enzyme_families_hit": len(unique_families),
            "hits_by_family": dict(sorted(families_hit.items(), key=lambda x: x[1], reverse=True)),
        },
        "interpretation_notes": [
            "unique_plastic_related_proteins_found counts unique reference protein accessions hit (best-hit homology).",
            "unique_reads_with_hits counts unique read IDs that produced at least one filtered hit.",
            "hits_by_family counts hit rows per enzyme family after filtering.",
            "estimated_reads_used is based on qc totals and run_params.max_reads if present; it is an estimate.",
        ],
    }

    write_json(args.out_json, summary)

    # Write Markdown
    ensure_dir(Path(args.out_md).parent)
    lines = []
    lines.append("# Run summary: data usage & plastic-enzyme hits\n\n")
    lines.append(f"- Generated (UTC): {summary['generated_utc']}\n\n")

    lines.append("## Data volume\n")
    lines.append(f"- Read1 total reads: {totals['read1_total_reads']:,}\n")
    lines.append(f"- Read2 total reads: {totals['read2_total_reads']:,}\n")
    lines.append(f"- Total reads (R1+R2): {totals['total_reads']:,}\n")
    lines.append(f"- Total bases (R1+R2): {totals['total_bases']:,}\n\n")

    lines.append("## Estimated data used in screening\n")
    lines.append(f"- max_reads param: {used['max_reads_param']}\n")
    lines.append(f"- Estimated used reads (R1): {used['estimated_used_reads_read1']:,}\n")
    lines.append(f"- Estimated used reads (R2): {used['estimated_used_reads_read2']:,}\n")
    lines.append(f"- Estimated used reads total: {used['estimated_used_reads_total']:,}\n")
    lines.append(f"- Note: {used['note']}\n\n")

    lines.append("## Plastic-related enzyme hits (homology-based)\n")
    lines.append(f"- Total hit rows: {total_hit_rows:,}\n")
    lines.append(f"- Reads with hits (unique read IDs): {len(unique_read_ids):,}\n")
    lines.append(f"- Plastic-related proteins found (unique accessions): {len(unique_protein_ids):,}\n")
    lines.append(f"- Enzyme families hit: {len(unique_families):,}\n\n")

    if families_hit:
        lines.append("### Hits by enzyme family\n")
        for fam, cnt in sorted(families_hit.items(), key=lambda x: x[1], reverse=True):
            lines.append(f"- {fam}: {cnt:,}\n")
        lines.append("\n")

    lines.append("## Output files\n")
    lines.append(f"- JSON: `{args.out_json}`\n")
    lines.append(f"- Markdown: `{args.out_md}`\n")

    Path(args.out_md).write_text("".join(lines), encoding="utf-8")

    print(f"[summary] Wrote: {args.out_json} and {args.out_md}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[summary] Error: {exc}", file=sys.stderr)
        sys.exit(1)
