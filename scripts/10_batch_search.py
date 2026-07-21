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


def prepare_subset(in_fastq, out_fastq, max_reads, trim_to=0):
    """Write the first ``max_reads`` records, optionally hard-trimmed to ``trim_to`` bp.

    Trimming exists to remove a batch confound: read length differs between the
    sequencing submissions in this BioProject (150 bp vs 159 bp). Because a blastx
    alignment cannot be longer than read_length/3 amino acids, a 150 bp read tops out
    at 50 aa while a 159 bp read reaches 53 aa. With a minimum-alignment-length filter
    anywhere near that ceiling, the shorter batch is systematically censored and looks
    depleted for reasons that have nothing to do with biology. Trimming every read to a
    common length gives all samples the same alignment-length ceiling.
    """
    ensure_dir(Path(out_fastq).parent)
    reads = 0
    with SEARCH.open_text(in_fastq) as fin, open(out_fastq, "w", encoding="utf-8") as fout:
        while reads < max_reads:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()
            if not qual:
                break
            if trim_to:
                seq = seq.rstrip("\n")[:trim_to] + "\n"
                qual = qual.rstrip("\n")[:trim_to] + "\n"
            fout.write(header)
            fout.write(seq)
            fout.write(plus)
            fout.write(qual)
            reads += 1
    return reads


def search_params(args):
    """Parameters that change the result; a cached sample is reused only if they match."""
    return {
        "depth_read_pairs": args.depth,
        "trim_to": args.trim_to,
        "min_aln_len": args.min_aln_len,
        "evalue": args.evalue,
        "pident": args.pident,
    }


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

    wanted = search_params(args)
    if out_path.exists() and out_path.stat().st_size > 0:
        # Reuse a cached result only if EVERY parameter that affects it matches.
        # Checking just the output's existence (or only the depth) would silently keep
        # a result produced under different settings, which is how a shallow smoke-test
        # run or an older filter threshold can contaminate a batch.
        previous = {}
        if meta_path.exists():
            try:
                previous = json.loads(meta_path.read_text(encoding="utf-8"))
            except (ValueError, OSError):
                previous = {}
        if all(previous.get(k) == v for k, v in wanted.items()):
            print(f"[skip] {run} already searched with these settings "
                  f"({count_hits(out_path):,} hits).")
            return out_path
        differing = [k for k, v in wanted.items() if previous.get(k) != v]
        print(f"[redo] {run} was searched with different settings "
              f"({', '.join(differing)}); re-running.")

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
            trim_note = f", trimmed to {args.trim_to} bp" if args.trim_to else ""
            print(f"  [{run}] subsampling {label} to {args.depth:,} reads{trim_note} ...",
                  flush=True)
            prepare_subset(str(fq), str(subset), args.depth, args.trim_to)

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
                    "reads_searched": args.depth * 2,
                    "hits": n_hits,
                    **search_params(args),
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
    parser.add_argument(
        "--trim-to", type=int, default=150,
        help="hard-trim every read to this many bp so all sequencing batches share the "
             "same alignment-length ceiling (0 disables trimming)",
    )
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--pident", type=float, default=30.0)
    parser.add_argument(
        "--min-aln-len", type=int, default=30,
        help="minimum alignment length in amino acids; must stay well below "
             "trim_to/3 or short-read batches get censored at the ceiling",
    )
    parser.add_argument("--limit", type=int, default=0, help="only process the first N runs (testing)")
    parser.add_argument("--commands-log", default="outputs/comparative/logs/commands.log")
    args = parser.parse_args()

    runs = read_metadata(args.metadata)
    if args.limit:
        runs = runs[: args.limit]

    # A blastx alignment cannot exceed read_length/3 amino acids. If the minimum
    # alignment length sits at or near that ceiling, hits pile up at the boundary and
    # samples with shorter reads look artificially depleted. Refuse to run in that regime.
    if args.trim_to:
        ceiling_aa = args.trim_to // 3
        if args.min_aln_len > ceiling_aa * 0.8:
            raise ValueError(
                f"--min-aln-len {args.min_aln_len} is too close to the {ceiling_aa} aa "
                f"ceiling implied by --trim-to {args.trim_to}. Alignments would be "
                f"censored at the boundary. Use --min-aln-len <= {int(ceiling_aa * 0.8)}."
            )

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
