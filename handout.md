# Handout — Windows → Linux port & research-output prep

> Living progress doc. If a session hits its limit, resume from **"Next steps"** below.
> Goal: make this metagenomics pipeline run on Linux (CachyOS), do everything up to
> the point of downloading the FASTQ dataset (owner will download that). Then produce
> research output (poster/conference paper) — later phase.

Last updated: 2026-07-18

---

## Project in one paragraph
Homology-based screen for plastic-degrading enzymes in a freshwater sediment shotgun
metagenome (NCBI SRA `SRR23872596`, ~27M read pairs, TezpurPTM_Brahmaputra). Pipeline:
metadata → QC → live-fetch UniProt enzyme refs → DIAMOND blastx homology search →
best-hit genus taxonomy → CPM aggregation → k-means clustering → network modularity
(Louvain) → markdown report. All Python. Reports *genetic potential only*, conservative
thesis-safe language. Key prior findings: AlkB most abundant (1,134 hits), *Aspergillus*
top genus (266 hits), 4 network modules.

## Repo shape
- `run_pipeline.py` — orchestrator, calls scripts 00→99 via subprocess.
- `scripts/00..99` — metadata, qc, fetch_refs, diamond_search, taxonomy, aggregate,
  cluster, network_modularity, run_summary, make_report. `utils.py` = shared helpers.
- `refs/` — cached UniProt enzyme FASTA + maps (already populated from prior run).
- `outputs/` — results/figures/tables from the **previous Windows run** (committed).
- Two PDFs + methodology.png + explanation.md = write-up assets.
- Remote: github.com/AbrarUI12/Metagenomic-Profiling-...  (branch: main)

## Environment findings (2026-07-18)
- OS: CachyOS (Arch-based), Python 3.14.6 system, **PEP-668 externally-managed** (no direct pip).
- Installed: `uv` (~/.local/bin/uv), numpy. Internet + 231G free disk OK.
- **Missing**: pandas, matplotlib, seaborn, scikit-learn, networkx, diamond, conda/mamba,
  blast, mmseqs, sra-tools (prefetch/fasterq-dump), pandoc.

## Windows-specific code to fix (identified)
1. `run_pipeline.py:16` — `Path("tools/diamond/diamond.exe")` (Windows binary path).
2. `scripts/03_diamond_search.py`:
   - `download_diamond()` (L86-104) downloads `diamond-windows.zip`.
   - `find_diamond()` (L107-114) only looks for `tools/diamond/diamond.exe`.
3. `.gitignore` ignores `*.exe` (harmless, but Linux diamond binary must not be ignored/committed).
4. `datetime.utcnow()` used throughout — deprecated in Py3.12+ (works, emits warning). Optional modernize.
5. `subprocess.run(shell=True)` — fine on Linux (uses /bin/sh, not fish).

## Decisions (confirmed by owner 2026-07-18)
- [x] Environment: **uv venv + static DIAMOND binary** (pinned Python 3.12).
- [x] FASTQ acquisition: **ENA direct download**, crash-safe — in-progress file restarts
      from scratch on power loss, finished files skipped (verified by size + md5).
- [x] Preserve prior Windows `outputs/` as baseline (kept at `outputs_windows_baseline/`, gitignored).

## Done (Linux port complete — everything up to dataset download)
- `.venv` created (Python 3.12) + all deps installed (`uv pip install -r requirements.txt`).
  pandas 3.0 / scikit-learn 1.9 / matplotlib 3.11 — all compatible.
- DIAMOND made **cross-platform**: `scripts/03_diamond_search.py` `download_diamond()` +
  `find_diamond()` now pick the right OS asset/binary name; `run_pipeline.py` diamond path too.
  DIAMOND v2.2.4 Linux binary downloaded to `tools/diamond/diamond` (gitignored via `tools/`).
- Modernized deprecated `datetime.utcnow()` → `datetime.now(timezone.utc)` in 6 files
  (utils, 00, 01, 02, 07_run_summary, 99). Output format unchanged.
- New `scripts/fetch_dataset.py`: crash-safe ENA downloader (+ gunzip). Logic + real
  interrupt/restart tested.
- README rewritten for Linux (uv venv, dataset download, DIAMOND auto/manual install).
- `.gitignore`: added `outputs/tmp/`, `outputs_windows_baseline/`.

### Verified on Linux (with .venv)
- 00 metadata (NCBI live) ✓ · 02 fetch_refs (UniProt live, 18,655 seqs) ✓
- 03 DIAMOND search ✓ (synthetic PETase read → 100% id hit, method=diamond)
- 06 cluster ✓ · 07 network ✓ (reproduced 4 modules, Q=0.255) · 99 report ✓
- Not separately run (low-risk, need real FASTQ / same libs already proven): 01 QC, 04 taxonomy, 05 aggregate.

### Note / gotcha
- `03_diamond_search.py` hardcodes intermediate paths `outputs/enzyme_hits/{raw,filtered}_{label}.tsv`
  (not overridable by `--out`). Running it clobbers those committed files. During testing this
  happened and was restored via `git checkout -- outputs/enzyme_hits/`. Harmless for the real
  single-dataset run (they get regenerated correctly), just can't isolate test runs there.

## Next steps (after this)
1. Owner commits the port (see commit commands provided in chat). Owner is author; assistant not co-author.
2. Owner downloads dataset: `.venv/bin/python scripts/fetch_dataset.py` (~1.7 GB gz → ~fastq).
3. Full run: `.venv/bin/python run_pipeline.py`. Compare against `outputs_windows_baseline/`.
4. Then: research-output phase (poster/conference paper) — separate effort.

## Commit protocol
- Assistant does NOT commit. Owner (Muhammad) commits; assistant only supplies terminal commands.
- Assistant is a contributor, not commit author.
