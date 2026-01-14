import argparse
import csv
import sys
import ssl
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts.utils import ensure_dir


def read_tsv(path):
    data = {}
    if not Path(path).exists():
        return data
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            data[row["protein_id"]] = row
    return data


def read_hit_proteins(path):
    proteins = set()
    if not Path(path).exists():
        return proteins
    with open(path, "r", encoding="utf-8") as f:
        next(f, None)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                proteins.add(parts[1])
    return proteins


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


def organism_to_genus(organism_name):
    if not organism_name:
        return ""
    return organism_name.split()[0]


def fetch_missing(accessions):
    rows = []
    if not accessions:
        return rows
    chunk = []
    for acc in accessions:
        chunk.append(acc)
        if len(chunk) == 100:
            rows.extend(fetch_chunk(chunk))
            chunk = []
    if chunk:
        rows.extend(fetch_chunk(chunk))
    return rows


def fetch_chunk(chunk):
    query = " OR ".join([f"accession:{acc}" for acc in chunk])
    tsv = uniprot_stream(
        query,
        "tsv",
        fields="accession,id,protein_name,organism_name,organism_id",
    )
    reader = csv.DictReader(tsv.splitlines(), delimiter="\t")
    rows = []
    for row in reader:
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--hits", default="outputs/enzyme_hits/enzyme_hits.tsv"
    )
    parser.add_argument(
        "--taxonomy", default="refs/protein_to_taxonomy.tsv"
    )
    args = parser.parse_args()

    taxonomy = read_tsv(args.taxonomy)
    hit_proteins = read_hit_proteins(args.hits)
    missing = [p for p in hit_proteins if p not in taxonomy]

    if missing:
        rows = fetch_missing(missing)
        for row in rows:
            acc = row.get("Entry") or row.get("accession")
            if not acc:
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

    ensure_dir(Path(args.taxonomy).parent)
    with open(args.taxonomy, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["protein_id", "taxid", "genus", "species", "source_db"])
        for acc, row in taxonomy.items():
            writer.writerow(
                [
                    row.get("protein_id", acc),
                    row.get("taxid", ""),
                    row.get("genus", ""),
                    row.get("species", ""),
                    row.get("source_db", "UniProt"),
                ]
            )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[taxonomy] Error: {exc}", file=sys.stderr)
        sys.exit(1)
