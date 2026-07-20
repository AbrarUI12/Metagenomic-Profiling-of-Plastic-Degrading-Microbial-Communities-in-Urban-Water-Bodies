"""Per-sample enzyme homology search across every run in the design table.

For each run listed in ``refs/sample_metadata.tsv`` this subsamples both mates to a
common depth, runs DIAMOND blastx against the same UniProt plastic-enzyme database
used by the single-sample pipeline, applies the same thresholds, and writes the
filtered hits to ``outputs/comparative/per_sample/<RUN>/enzyme_hits.tsv``.

Equal-depth subsampling matters: library sizes here range 13.5 M-70 M read pairs, so
raw hit counts are otherwise dominated by how deeply a sample happened to be
sequenced rather than by its biology.

**Resumable.** A run whose output already exists (and is non-empty) is skipped, so an
interrupted batch is continued by re-issuing the same command. Output is written to a
``.part`` file and renamed only on success, so a killed process never leaves a
half-written table that would later be mistaken for a finished one.

This never writes to the single-sample paths under ``outputs/enzyme_hits/`` — those
are hardcoded in ``03_diamond_search.py`` and would be clobbered.
"""

import argparse
import csv
import importlib.util
import json
import shutil
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts.utils import ensure_dir, log_command, run_cmd


def _load_search_module():
    """Import 03_diamond_search.py, whose module name cannot be imported normally."""
    path = Path(__file__).resolve().parent / "03_diamond_search.py"
    spec = importlib.util.spec_from_file_location("diamond_search", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SEARCH = _load_search_module()


def read_metadata(path):
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def build_diamond_db(diamond, fasta_db, commands_log):
    """Build the DIAMOND index once per batch (seconds for ~18k proteins).

    Done here rather than per sample so 18 runs don't rebuild the same index 36
    times, and so the index is guaranteed to match the binary in use.
    """
    db_prefix = str(Path(fasta_db).with_suffix(""))
    cmd = f"\"{diamond}\" makedb --in \"{fasta_db}\" -d \"{db_prefix}\""
    run_cmd(cmd, commands_log=commands_log)
    return db_prefix


def diamond_blastx(diamond, db_prefix, fastq, out_tsv, threads, evalue, pident, commands_log):
    cmd = (
        f"\"{diamond}\" blastx -d \"{db_prefix}\" -q \"{fastq}\" "
        f"-o \"{out_tsv}\" --evalue {evalue} --id {pident} "
        f"--max-target-seqs 1 --threads {threads} "
        "--outfmt 6 qseqid sseqid evalue bitscore pident length"
    )
    run_cmd(cmd, commands_log=commands_log)


def count_hits(path):
    with open(path, "r", encoding="utf-8") as f:
        next(f, None)
        return sum(1 for line in f if line.strip())


def process_sample(run, args, diamond, db_prefix, family_map):
    """Search one run. Returns the path to its per-sample hits table."""
    out_dir = Path(args.out_dir) / run
    out_path = out_dir / "enzyme_hits.tsv"
    meta_path = out_dir / "sample_meta.json"

    if out_path.exists() and out_path.stat().st_size > 0:
        # Only reuse a result that was produced at the SAME depth. Otherwise a
        # leftover from a shallow smoke-test run would be silently kept in a
        # full-depth batch, reintroducing exactly the unequal-depth artifact that
        # equal-depth subsampling exists to remove.
        previous_depth = None
        if meta_path.exists():
            try:
                previous_depth = json.loads(meta_path.read_text(encoding="utf-8")).get(
                    "depth_read_pairs"
                )
            except (ValueError, OSError):
                previous_depth = None
        if previous_depth == args.depth:
            print(f"[skip] {run} already searched at this depth ({count_hits(out_path):,} hits).")
            return out_path
        print(
            f"[redo] {run} was searched at depth {previous_depth} "
            f"(now {args.depth:,}); re-running."
        )

    fq1 = Path(args.data_dir) / f"{run}_1.fastq.gz"
    fq2 = Path(args.data_dir) / f"{run}_2.fastq.gz"
    for fq in (fq1, fq2):
        if not fq.exists():
            raise FileNotFoundError(f"{fq} missing - run scripts/fetch_dataset.py first")

    ensure_dir(out_dir)
    tmp_dir = Path(args.tmp_dir) / run
    ensure_dir(tmp_dir)

    filtered_paths = []
    try:
        for label, fq in [("read1", fq1), ("read2", fq2)]:
            subset = tmp_dir / f"{label}_subset.fastq"
            print(f"  [{run}] subsampling {label} to {args.depth:,} reads ...", flush=True)
            SEARCH.prepare_fastq_subset(str(fq), str(subset), args.depth)

            raw_out = tmp_dir / f"raw_{label}.tsv"
            filtered_out = tmp_dir / f"filtered_{label}.tsv"
            print(f"  [{run}] DIAMOND blastx {label} ...", flush=True)
            diamond_blastx(
                diamond, db_prefix, str(subset), str(raw_out),
                args.threads, args.evalue, args.pident, args.commands_log,
            )
            SEARCH.filter_hits(
                str(raw_out), str(filtered_out), family_map,
                args.evalue, args.pident, args.min_aln_len,
            )
            filtered_paths.append(str(filtered_out))
            # Free the subset immediately; at 13 M reads each is several GB.
            subset.unlink(missing_ok=True)
            raw_out.unlink(missing_ok=True)

        part_path = out_path.with_suffix(".tsv.part")
        SEARCH.merge_outputs(filtered_paths, str(part_path))
        part_path.rename(out_path)

        n_hits = count_hits(out_path)
        meta_path.write_text(
            json.dumps(
                {
                    "run_accession": run,
                    "depth_read_pairs": args.depth,
                    "reads_searched": args.depth * 2,
                    "hits": n_hits,
                    "evalue": args.evalue,
                    "pident": args.pident,
                    "min_aln_len": args.min_aln_len,
                },
                indent=2,
                sort_keys=True,
            ),
            encoding="utf-8",
        )
        print(f"[done] {run}: {n_hits:,} filtered hits.")
        return out_path
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default="refs/sample_metadata.tsv")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument("--out-dir", default="outputs/comparative/per_sample")
    parser.add_argument("--tmp-dir", default="outputs/comparative/tmp")
    parser.add_argument("--fasta-db", default="refs/plastic_enzymes.fasta")
    parser.add_argument("--family-map", default="refs/enzyme_family_map.tsv")
    parser.add_argument(
        "--depth", type=int, default=13_000_000,
        help="read pairs kept per sample (default 13M ~ the smallest library)",
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--pident", type=float, default=30.0)
    parser.add_argument("--min-aln-len", type=int, default=50)
    parser.add_argument("--limit", type=int, default=0, help="only process the first N runs (testing)")
    parser.add_argument("--commands-log", default="outputs/comparative/logs/commands.log")
    args = parser.parse_args()

    runs = read_metadata(args.metadata)
    if args.limit:
        runs = runs[: args.limit]

    diamond = SEARCH.find_diamond()
    if not diamond:
        raise RuntimeError(
            "DIAMOND not found. Install it or let run_pipeline.py download it into tools/diamond/."
        )
    family_map = SEARCH.read_family_map(args.family_map)

    ensure_dir(Path(args.commands_log).parent)
    log_command(args.commands_log, f"BATCH_SEARCH depth={args.depth} samples={len(runs)}")
    db_prefix = build_diamond_db(diamond, args.fasta_db, args.commands_log)

    print(f"Searching {len(runs)} sample(s) at {args.depth:,} read pairs each.\n")
    for idx, row in enumerate(runs, 1):
        run = row["run_accession"]
        print(f"[{idx}/{len(runs)}] {run}  ({row['site']}, {row['season']})")
        process_sample(run, args, diamond, db_prefix, family_map)

    print(f"\nPer-sample hits written under {args.out_dir}/")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[batch_search] Error: {exc}", file=sys.stderr)
        sys.exit(1)
