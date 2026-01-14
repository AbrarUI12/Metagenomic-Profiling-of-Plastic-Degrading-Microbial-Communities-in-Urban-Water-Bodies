import argparse
import csv
import json
import sys
import ssl
from datetime import datetime
from pathlib import Path
from urllib.request import urlopen

sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts.utils import ensure_dir, write_json


def read_csv(path):
    with open(path, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def fetch_runinfo_csv(srr):
    esearch = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=sra&term={srr}&retmode=json"
    )
    try:
        with urlopen(esearch) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except Exception:
        ctx = ssl._create_unverified_context()
        with urlopen(esearch, context=ctx) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    ids = data.get("esearchresult", {}).get("idlist", [])
    if not ids:
        raise RuntimeError("No SRA IDs found for accession")
    sra_id = ids[0]
    efetch = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=sra&id={sra_id}&rettype=runinfo&retmode=text"
    )
    try:
        with urlopen(efetch) as resp:
            return resp.read().decode("utf-8")
    except Exception:
        ctx = ssl._create_unverified_context()
        with urlopen(efetch, context=ctx) as resp:
            return resp.read().decode("utf-8")


def choose_row(rows, srr):
    if not rows:
        return None
    for row in rows:
        if row.get("Run") == srr:
            return row
    return rows[0]


def summarize_metadata(row):
    if not row:
        return {}
    fields = [
        "Run",
        "BioSample",
        "BioProject",
        "SampleName",
        "SampleTitle",
        "Sample",
        "ScientificName",
        "LibraryStrategy",
        "LibrarySource",
        "LibrarySelection",
        "Platform",
        "Instrument",
        "ReleaseDate",
        "LoadDate",
        "spots",
        "bases",
        "avgLength",
    ]
    summary = {}
    for f in fields:
        if f in row and row[f]:
            summary[f] = row[f]
    return summary


def write_markdown(out_md, summary, source_note, out_json_path):
    lines = []
    lines.append("# Metadata summary\n")
    lines.append(f"- Retrieval date (UTC): {datetime.utcnow().isoformat()}Z\n")
    lines.append(f"- Source: {source_note}\n")
    if summary:
        lines.append("\n## Key fields\n")
        for k, v in summary.items():
            lines.append(f"- {k}: {v}\n")
    else:
        lines.append("\nNo metadata fields were found.\n")
    lines.append("\n## Notes\n")
    lines.append(
        "- This summary is derived from SRA run metadata and may not include "
        "detailed environmental context beyond what submitters provide.\n"
    )
    lines.append(f"- Full metadata JSON: `{out_json_path}`\n")

    ensure_dir(Path(out_md).parent)
    with open(out_md, "w", encoding="utf-8") as f:
        f.writelines(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--srr", default="SRR23872596")
    parser.add_argument("--sra-run-table", default=None)
    parser.add_argument("--out-md", default="outputs/metadata_summary.md")
    parser.add_argument("--out-json", default="outputs/metadata/metadata.json")
    parser.add_argument("--out-runinfo", default="outputs/metadata/sra_runinfo.csv")
    parser.add_argument("--commands-log", default=None)
    args = parser.parse_args()

    rows = []
    source_note = ""
    if args.sra_run_table and Path(args.sra_run_table).exists():
        rows = read_csv(args.sra_run_table)
        source_note = f"Local SraRunTable.csv ({args.sra_run_table})"
    else:
        runinfo_csv = fetch_runinfo_csv(args.srr)
        ensure_dir(Path(args.out_runinfo).parent)
        Path(args.out_runinfo).write_text(runinfo_csv, encoding="utf-8")
        rows = list(csv.DictReader(runinfo_csv.splitlines()))
        source_note = f"NCBI SRA runinfo (SRR={args.srr})"

    row = choose_row(rows, args.srr)
    summary = summarize_metadata(row)
    metadata = {
        "srr": args.srr,
        "source": source_note,
        "row": row,
        "summary": summary,
    }
    write_json(args.out_json, metadata)
    write_markdown(args.out_md, summary, source_note, args.out_json)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[metadata] Error: {exc}", file=sys.stderr)
        sys.exit(1)
