"""Orchestrator for the spatio-temporal comparative analysis (steps 10 -> 14).

Runs every sample in ``refs/sample_metadata.tsv`` through the enzyme homology search
at a common depth, then builds cross-sample matrices, comparative statistics, figures
and a markdown report.

This is separate from ``run_pipeline.py`` (the original single-sample analysis), which
it neither calls nor modifies. All output lands under ``outputs/comparative/``.

Usage
-----
    # Full analysis (hours: 18 samples x DIAMOND blastx)
    .venv/bin/python run_comparative.py --threads 12

    # Quick end-to-end smoke test (tiny depth, 2 samples)
    .venv/bin/python run_comparative.py --threads 12 --depth 200000 --limit 2

The search step is resumable: samples already searched are skipped, so re-issuing the
same command after an interruption continues where it stopped.
"""

import argparse
import subprocess
import sys
from pathlib import Path

from scripts.utils import ensure_dir, log_command

STEPS_DIR = Path("scripts")
COMMANDS_LOG = "outputs/comparative/logs/commands.log"


def run_step(python, script, step_args, commands_log):
    cmd = [python, str(STEPS_DIR / script)] + step_args
    printable = " ".join(cmd)
    print(f"\n=== {script} ===", flush=True)
    log_command(commands_log, printable)
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise RuntimeError(f"{script} failed with exit code {result.returncode}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default="refs/sample_metadata.tsv")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument("--out-root", default="outputs/comparative")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--depth", type=int, default=13_000_000,
        help="read pairs searched per sample (default 13M ~ smallest library)",
    )
    parser.add_argument(
        "--trim-to", type=int, default=150,
        help="hard-trim reads to this length so all sequencing batches are comparable",
    )
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--pident", type=float, default=30.0)
    parser.add_argument("--min-aln-len", type=int, default=30)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--limit", type=int, default=0, help="only process first N samples (testing)")
    parser.add_argument("--skip-search", action="store_true",
                        help="reuse existing per-sample hits; only rebuild matrices/stats/figures")
    args = parser.parse_args()

    out_root = Path(args.out_root)
    tables = out_root / "tables"
    figures = out_root / "figures"
    per_sample = out_root / "per_sample"
    for d in [out_root, tables, figures, per_sample, out_root / "logs", out_root / "tmp"]:
        ensure_dir(d)

    if not Path(args.metadata).exists():
        raise FileNotFoundError(
            f"{args.metadata} not found - run scripts/make_sample_metadata.py first."
        )

    python = sys.executable
    log_command(COMMANDS_LOG, f"RUN_COMPARATIVE depth={args.depth} threads={args.threads}")

    if not args.skip_search:
        search_args = [
            "--metadata", args.metadata,
            "--data-dir", args.data_dir,
            "--out-dir", str(per_sample),
            "--tmp-dir", str(out_root / "tmp"),
            "--depth", str(args.depth),
            "--threads", str(args.threads),
            "--trim-to", str(args.trim_to),
            "--evalue", str(args.evalue),
            "--pident", str(args.pident),
            "--min-aln-len", str(args.min_aln_len),
            "--commands-log", COMMANDS_LOG,
        ]
        if args.limit:
            search_args += ["--limit", str(args.limit)]
        run_step(python, "10_batch_search.py", search_args, COMMANDS_LOG)

    matrix_args = [
        "--metadata", args.metadata,
        "--per-sample-dir", str(per_sample),
        "--out-dir", str(tables),
        "--depth", str(args.depth),
    ]
    if args.limit:
        matrix_args += ["--limit", str(args.limit)]
    run_step(python, "11_build_matrices.py", matrix_args, COMMANDS_LOG)

    run_step(python, "12_comparative_stats.py", [
        "--tables-dir", str(tables),
        "--out-dir", str(tables),
        "--permutations", str(args.permutations),
        "--seed", str(args.seed),
    ], COMMANDS_LOG)

    run_step(python, "13_comparative_plots.py", [
        "--tables-dir", str(tables),
        "--fig-dir", str(figures),
    ], COMMANDS_LOG)

    run_step(python, "14_comparative_report.py", [
        "--tables-dir", str(tables),
        "--fig-dir", str(figures),
        "--out-md", str(out_root / "report.md"),
    ], COMMANDS_LOG)

    print(f"\nComparative analysis complete. Report: {out_root / 'report.md'}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[interrupted] Re-run the same command; searched samples are skipped.")
        sys.exit(130)
    except Exception as exc:
        print(f"[run_comparative] Error: {exc}", file=sys.stderr)
        sys.exit(1)
