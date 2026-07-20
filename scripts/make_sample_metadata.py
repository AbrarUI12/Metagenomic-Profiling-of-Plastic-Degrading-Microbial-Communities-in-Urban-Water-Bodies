"""Build the committed sample design table for the comparative analysis.

Pulls run metadata for the BioProject from ENA and derives the two experimental
factors the comparative analysis tests: **site** (spatial position along the
Brahmaputra) and **season** (sampling campaign). The result is written to
``refs/sample_metadata.tsv`` and committed, so the analysis is deterministic and
does not depend on a live ENA call at analysis time.

Re-run this only to refresh the design table:

    .venv/bin/python scripts/make_sample_metadata.py
"""

import argparse
import csv
import sys
from pathlib import Path
from urllib.request import Request, urlopen

sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts.utils import ensure_dir

ENA_FILEREPORT = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession={acc}&result=read_run"
    "&fields=run_accession,sample_title,sample_alias,location,lat,lon,country,"
    "collection_date,library_strategy,base_count,read_count,scientific_name"
    "&format=tsv"
)

# Uzanbazar is the river frontage of Guwahati (~1.1 M people), by far the largest
# urban centre on the Brahmaputra in this transect. Every other site is a town or
# rural reach, so "urban" is a single-site contrast and is reported as such.
URBAN_SITES = {"Uzanbazar"}


def season_from_date(collection_date):
    """Map a collection date to the Assam hydrological season.

    Sampling campaigns in this BioProject fall in Sep 2021, Dec 2021-Jan 2022 and
    Apr 2022. Seasons are derived from the date (ground truth) rather than the
    sample-code suffixes, which are ambiguous.
    """
    if not collection_date:
        return "unknown"
    try:
        month = int(str(collection_date).split("-")[1])
    except (IndexError, ValueError):
        return "unknown"
    if month in (6, 7, 8):
        return "monsoon"
    if month in (9, 10, 11):
        return "post_monsoon"
    if month in (12, 1, 2):
        return "winter"
    return "pre_monsoon"


def site_from_country(country, sample_title):
    """ENA's `country` field looks like 'India: Tezpur' -> site name 'Tezpur'."""
    if country and ":" in country:
        return country.split(":", 1)[1].strip()
    if country:
        return country.strip()
    return (sample_title or "unknown").strip()


def fetch_runs(accession):
    url = ENA_FILEREPORT.format(acc=accession)
    req = Request(url, headers={"User-Agent": "plastic-metagenome-pipeline/1.0"})
    with urlopen(req, timeout=60) as resp:
        text = resp.read().decode("utf-8")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        raise RuntimeError(f"ENA returned no run records for {accession}")
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]


def build_rows(runs):
    rows = []
    for r in runs:
        site = site_from_country(r.get("country", ""), r.get("sample_title", ""))
        collection_date = r.get("collection_date", "")
        read_count = r.get("read_count", "") or "0"
        rows.append(
            {
                "run_accession": r.get("run_accession", ""),
                "site": site,
                "site_code": r.get("sample_title", ""),
                "lat": r.get("lat", ""),
                "lon": r.get("lon", ""),
                "season": season_from_date(collection_date),
                "collection_date": collection_date,
                "urban": 1 if site in URBAN_SITES else 0,
                "read_pairs_total": int(read_count),
                "sample_alias": r.get("sample_alias", ""),
            }
        )
    # Sort upstream -> downstream (descending latitude), then by date, so the table
    # reads along the river the same way the figures do.
    rows.sort(key=lambda x: (-float(x["lat"] or 0), x["collection_date"]))
    return rows


FIELDS = [
    "run_accession",
    "site",
    "site_code",
    "lat",
    "lon",
    "season",
    "collection_date",
    "urban",
    "read_pairs_total",
    "sample_alias",
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acc", default="PRJNA944918", help="BioProject accession")
    parser.add_argument("--out", default="refs/sample_metadata.tsv")
    args = parser.parse_args()

    print(f"Querying ENA for {args.acc} ...")
    rows = build_rows(fetch_runs(args.acc))

    ensure_dir(Path(args.out).parent)
    with open(args.out, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    sites = sorted({r["site"] for r in rows})
    seasons = sorted({r["season"] for r in rows})
    print(f"Wrote {len(rows)} runs -> {args.out}")
    print(f"  sites   ({len(sites)}): {', '.join(sites)}")
    print(f"  seasons ({len(seasons)}): {', '.join(seasons)}")
    for season in seasons:
        n = sum(1 for r in rows if r["season"] == season)
        print(f"    {season:14} {n} samples")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[sample_metadata] Error: {exc}", file=sys.stderr)
        sys.exit(1)
