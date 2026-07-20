"""Combine the per-sample hit tables into cross-sample matrices.

Turns the 18 independent ``enzyme_hits.tsv`` files produced by ``10_batch_search.py``
into the sample-by-feature matrices the comparative statistics operate on:

* ``sample_enzyme_matrix_{raw,cpm}.csv`` - rows = samples, cols = enzyme families
* ``sample_genus_matrix_{raw,cpm}.csv``  - rows = samples, cols = genera
* ``sample_summary.tsv``                 - per-sample totals, diversity and design factors

Best-hit taxonomy is resolved once for the union of proteins hit across all samples
(one UniProt round instead of 18), reusing the cache in ``refs/protein_to_taxonomy.tsv``.

CPM is counts per million *searched reads*. Because every sample was searched at the
same depth, the raw and CPM matrices are proportional here; CPM is kept so the numbers
stay comparable to the single-sample report and remain correct if depth ever varies.
"""

import argparse
import csv
import importlib.util
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from scripts.utils import ensure_dir


def _load_taxonomy_module():
    path = Path(__file__).resolve().parent / "04_taxonomy_map.py"
    spec = importlib.util.spec_from_file_location("taxonomy_map", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


TAX = _load_taxonomy_module()


def read_metadata(path):
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_hits(per_sample_dir, runs):
    """Load every per-sample hit table into one long dataframe with a run column."""
    frames = []
    missing = []
    for run in runs:
        path = Path(per_sample_dir) / run / "enzyme_hits.tsv"
        if not path.exists():
            missing.append(run)
            continue
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            # A sample can legitimately yield zero hits; keep it as an all-zero row.
            df = pd.DataFrame(columns=["read_id", "protein_id", "family"])
        df["run_accession"] = run
        frames.append(df)
    if missing:
        raise FileNotFoundError(
            f"No hit table for {len(missing)} run(s): {', '.join(missing)}. "
            "Run scripts/10_batch_search.py first."
        )
    return pd.concat(frames, ignore_index=True)


def resolve_taxonomy(hits, taxonomy_path):
    """Ensure every hit protein has a genus, fetching only what the cache lacks."""
    taxonomy = TAX.read_tsv(taxonomy_path)
    hit_proteins = set(hits["protein_id"].dropna().astype(str))
    missing = sorted(p for p in hit_proteins if p not in taxonomy)

    if missing:
        print(f"  resolving {len(missing)} new protein->genus mappings from UniProt ...")
        for row in TAX.fetch_missing(missing):
            acc = row.get("Entry") or row.get("accession")
            if not acc:
                continue
            org = row.get("Organism") or row.get("organism_name") or ""
            taxonomy[acc] = {
                "protein_id": acc,
                "taxid": row.get("Organism ID") or row.get("organism_id") or "",
                "genus": TAX.organism_to_genus(org),
                "species": org,
                "source_db": "UniProt",
            }
        ensure_dir(Path(taxonomy_path).parent)
        with open(taxonomy_path, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["protein_id", "taxid", "genus", "species", "source_db"])
            for acc, row in taxonomy.items():
                writer.writerow([
                    row.get("protein_id", acc), row.get("taxid", ""),
                    row.get("genus", ""), row.get("species", ""),
                    row.get("source_db", "UniProt"),
                ])

    genus_of = {acc: (row.get("genus") or "Unassigned") for acc, row in taxonomy.items()}
    hits["genus"] = hits["protein_id"].map(genus_of).fillna("Unassigned")
    return hits


def pivot_counts(hits, column, runs):
    """samples x category count matrix, with every run present even if it had no hits."""
    if hits.empty or column not in hits:
        return pd.DataFrame(index=pd.Index(runs, name="run_accession"))
    matrix = (
        hits.groupby(["run_accession", column])
        .size()
        .unstack(fill_value=0)
    )
    matrix = matrix.reindex(runs, fill_value=0)
    matrix.index.name = "run_accession"
    return matrix.astype(int)


def shannon(counts):
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total <= 0:
        return 0.0
    p = counts[counts > 0] / total
    return float(-(p * np.log(p)).sum())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default="refs/sample_metadata.tsv")
    parser.add_argument("--per-sample-dir", default="outputs/comparative/per_sample")
    parser.add_argument("--taxonomy", default="refs/protein_to_taxonomy.tsv")
    parser.add_argument("--out-dir", default="outputs/comparative/tables")
    parser.add_argument(
        "--depth", type=int, default=13_000_000,
        help="read pairs searched per sample (must match 10_batch_search.py)",
    )
    parser.add_argument(
        "--limit", type=int, default=0,
        help="only use the first N runs (must match 10_batch_search.py when testing)",
    )
    args = parser.parse_args()

    meta = pd.DataFrame(read_metadata(args.metadata))
    if args.limit:
        meta = meta.head(args.limit)
    runs = meta["run_accession"].tolist()

    print(f"Loading hits for {len(runs)} samples ...")
    hits = load_hits(args.per_sample_dir, runs)
    hits = resolve_taxonomy(hits, args.taxonomy)

    out_dir = Path(args.out_dir)
    ensure_dir(out_dir)

    reads_searched = args.depth * 2  # both mates were searched
    results = {}
    for name, column in [("enzyme", "family"), ("genus", "genus")]:
        raw = pivot_counts(hits, column, runs)
        cpm = raw / reads_searched * 1e6
        raw.to_csv(out_dir / f"sample_{name}_matrix_raw.csv")
        cpm.to_csv(out_dir / f"sample_{name}_matrix_cpm.csv")
        results[name] = raw
        print(f"  sample x {name}: {raw.shape[0]} samples x {raw.shape[1]} features")

    enzyme_raw, genus_raw = results["enzyme"], results["genus"]
    summary = meta.set_index("run_accession").copy()
    summary["total_hits"] = enzyme_raw.sum(axis=1)
    summary["total_hits_cpm"] = summary["total_hits"] / reads_searched * 1e6
    summary["genus_richness"] = (genus_raw > 0).sum(axis=1)
    summary["genus_shannon"] = [shannon(genus_raw.loc[r]) for r in summary.index]
    summary["family_richness"] = (enzyme_raw > 0).sum(axis=1)
    summary["family_shannon"] = [shannon(enzyme_raw.loc[r]) for r in summary.index]
    for family in enzyme_raw.columns:
        summary[f"hits_{family}"] = enzyme_raw[family]
    summary["reads_searched"] = reads_searched
    summary.to_csv(out_dir / "sample_summary.tsv", sep="\t")

    print(f"\nTotal hits across all samples: {int(summary['total_hits'].sum()):,}")
    print(f"Wrote matrices and sample_summary.tsv to {out_dir}/")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[build_matrices] Error: {exc}", file=sys.stderr)
        sys.exit(1)
