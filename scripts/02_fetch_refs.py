import argparse
import csv
import sys
import ssl
from datetime import datetime
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts.utils import ensure_dir


FAMILIES = {
    "PETase_like": {
        "query": '(protein_name:"poly(ethylene terephthalate) hydrolase" OR protein_name:PETase)',
        "description": "PETase/cutinase-like PET hydrolases",
    },
    "MHETase_like": {
        "query": '(protein_name:MHETase OR protein_name:"mono(2-hydroxyethyl) terephthalate hydrolase")',
        "description": "MHETase-like hydrolases",
    },
    "Cutinase_like": {
        "query": "protein_name:cutinase",
        "description": "Cutinase-like esterases",
    },
    "AlkB": {
        "query": '(protein_name:"alkane 1-monooxygenase" OR gene:alkB)',
        "description": "Alkane monooxygenase (AlkB) family",
    },
    "Polyesterase_like": {
        "query": '(protein_name:"polyesterase" OR protein_name:"polyester hydrolase" OR protein_name:"polyurethane esterase")',
        "description": "Polyesterase/polyurethane esterase-like enzymes",
    },
}


def uniprot_stream(query, fmt, fields=None):
    base = "https://rest.uniprot.org/uniprotkb/stream"
    params = {"query": query, "format": fmt}
    if fields:
        params["fields"] = fields
    url = base + "?" + urlencode(params)
    req = Request(url, headers={"User-Agent": "plastic-metagenome-pipeline/1.0"})
    try:
        with urlopen(req) as resp:
            return resp.read().decode("utf-8")
    except Exception:
        ctx = ssl._create_unverified_context()
        with urlopen(req, context=ctx) as resp:
            return resp.read().decode("utf-8")


def parse_fasta(fasta_text):
    seqs = {}
    header = None
    chunks = []
    for line in fasta_text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                seqs[header] = "".join(chunks)
            header = line[1:]
            chunks = []
        else:
            chunks.append(line)
    if header:
        seqs[header] = "".join(chunks)
    return seqs


def extract_accession(header):
    if "|" in header:
        parts = header.split("|")
        if len(parts) >= 2:
            return parts[1]
    return header.split()[0]


def fasta_records_with_accession(fasta_text):
    seqs = parse_fasta(fasta_text)
    records = {}
    for header, seq in seqs.items():
        acc = extract_accession(header)
        if acc not in records:
            records[acc] = (header, seq)
    return records


def parse_tsv(tsv_text):
    rows = []
    reader = csv.DictReader(tsv_text.splitlines(), delimiter="\t")
    for row in reader:
        rows.append(row)
    return rows


def organism_to_genus(organism_name):
    if not organism_name:
        return ""
    return organism_name.split()[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-fasta", default="refs/plastic_enzymes.fasta")
    parser.add_argument("--out-map", default="refs/enzyme_family_map.tsv")
    parser.add_argument("--out-taxonomy", default="refs/protein_to_taxonomy.tsv")
    parser.add_argument("--out-readme", default="refs/README_refs.md")
    parser.add_argument("--raw-dir", default="refs/raw")
    args = parser.parse_args()

    ensure_dir(args.raw_dir)

    all_records = {}
    family_map = {}
    taxonomy = {}
    family_counts = {}

    for family, meta in FAMILIES.items():
        query = meta["query"]
        fasta_text = uniprot_stream(query, "fasta")
        tsv_text = uniprot_stream(
            query,
            "tsv",
            fields="accession,id,protein_name,organism_name,organism_id",
        )

        raw_fasta_path = Path(args.raw_dir) / f"uniprot_{family}.fasta"
        raw_tsv_path = Path(args.raw_dir) / f"uniprot_{family}.tsv"
        raw_fasta_path.write_text(fasta_text, encoding="utf-8")
        raw_tsv_path.write_text(tsv_text, encoding="utf-8")

        records = fasta_records_with_accession(fasta_text)
        family_counts[family] = len(records)

        for acc, (header, seq) in records.items():
            if acc not in all_records:
                all_records[acc] = (header, seq)
            if acc not in family_map:
                family_map[acc] = family

        tsv_rows = parse_tsv(tsv_text)
        for row in tsv_rows:
            acc = row.get("Entry") or row.get("accession")
            if not acc:
                continue
            if acc in taxonomy:
                continue
            org = row.get("Organism") or row.get("organism_name") or ""
            taxid = row.get("Organism ID") or row.get("organism_id") or ""
            taxonomy[acc] = {
                "protein_id": acc,
                "taxid": taxid,
                "genus": organism_to_genus(org),
                "species": org,
                "source_db": "UniProt",
            }

    ensure_dir(Path(args.out_fasta).parent)
    with open(args.out_fasta, "w", encoding="utf-8") as f:
        for acc, (header, seq) in all_records.items():
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")

    with open(args.out_map, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["protein_id", "family"])
        for acc, family in family_map.items():
            writer.writerow([acc, family])

    with open(args.out_taxonomy, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["protein_id", "taxid", "genus", "species", "source_db"])
        for acc, row in taxonomy.items():
            writer.writerow(
                [
                    row["protein_id"],
                    row["taxid"],
                    row["genus"],
                    row["species"],
                    row["source_db"],
                ]
            )

    retrieval_date = datetime.utcnow().strftime("%Y-%m-%d")
    lines = []
    lines.append("# Reference sequences for plastic-degrading enzymes\n")
    lines.append(f"- Retrieval date (UTC): {retrieval_date}\n")
    lines.append("- Source: UniProtKB REST API (stream endpoint)\n")
    lines.append("\n## Queries\n")
    for family, meta in FAMILIES.items():
        lines.append(f"- {family}: {meta['query']}\n")
    lines.append("\n## Counts\n")
    for family, count in family_counts.items():
        lines.append(f"- {family}: {count} sequences\n")
    lines.append("\n## Files\n")
    lines.append(f"- Combined FASTA: `{args.out_fasta}`\n")
    lines.append(f"- Family map: `{args.out_map}`\n")
    lines.append(f"- Taxonomy map: `{args.out_taxonomy}`\n")
    lines.append(f"- Raw downloads: `{args.raw_dir}`\n")

    Path(args.out_readme).write_text("".join(lines), encoding="utf-8")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[refs] Error: {exc}", file=sys.stderr)
        sys.exit(1)
