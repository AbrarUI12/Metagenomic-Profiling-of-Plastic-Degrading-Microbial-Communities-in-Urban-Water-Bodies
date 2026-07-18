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

## Full dataset inventory (BioProject PRJNA944918 — Brahmaputra river sediment metagenome)
Confirmed via ENA: **18 runs**, all "sediment metagenome" from Brahmaputra locations
(Tezpur, Morigaon, Dhubri, Dibrugarh, Uzanbazar, Sadiya, Palasbari, Dhuburi, Tsk, Sad...).
Original work used only **SRR23872596** (Tezpur, TZB).

- Total **compressed download (.fastq.gz): ~55.5 GB**
- Total **uncompressed (.fastq): ~260 GB** (estimate; could be ~285 GB)
- Free disk on /home: **~231 GB** → **all-uncompressed-at-once does NOT fit.**
- Largest single run: SRR28925922 (7.6 GB gz → ~42 GB fastq).

### AMPLICON-vs-shotgun — RESOLVED (2026-07-19)
SRA tags every run `library_strategy = AMPLICON` (selection PCR), which conflicted with
`explanation.md` ("shotgun metagenomics"). Resolved definitively:
- **Peer-reviewed publication** for PRJNA944918: *"Metagenomic insights into microbial
  community, functional annotation, and antibiotic resistance genes in Himalayan Brahmaputra
  River sediment, India"*, Frontiers in Microbiology 2024 (PMC11614985) — states
  **"high-throughput shotgun metagenomics."**
- **Empirical:** baseline SRR23872596 hits span 4 unrelated enzyme families at high identity
  (469 hits ≥70% pident; median 55.8%) — impossible from 16S amplicon reads. → data is shotgun.
- **Conclusion:** the SRA `AMPLICON` label is a **mislabel**. For the paper, describe as shotgun
  metagenomics, cite PMC11614985, and do NOT copy the SRA strategy tag.
- Note: the published paper covers **6 samples** (BRS 1-6 = the PTM series: SRR23872593-597 +
  SRR23961995). The other 12 runs were added to the same BioProject later.

### Storage / HDD (confirmed 2026-07-19)
- HDD `/mnt/shinrinyoku`: **NTFS (ntfs-3g/fuseblk), ~455 GB free**, writable as user. Fits all.
- Plan: DATA on HDD; keep `.venv` + `tools/diamond` on the SSD (/home). Peak (gz+fastq) ~315 GB — fits.
- `scripts/fetch_dataset.py` now accepts `--acc PRJNA944918` (whole project) or a single run;
  same crash-safe per-file logic. Verified: lists 18 runs / 36 files / 55.5 GB and writes to HDD.

**Processing approach — DECIDED: gzip-native (2026-07-19).** The pipeline now reads `.fastq.gz`
directly, so no uncompression is ever needed — total footprint is just the ~55 GB of .gz.
- Added `open_text()` to `utils.py` (transparent gz/plain reader); `01_qc_sanity.py`,
  `03` subset + kmer-fallback now use it. DIAMOND reads .gz natively. Verified on synthetic gz.
- `.gitignore` now ignores `*.fastq.gz` too.
- Because it's only 55 GB and gz reads are lossless-identical, **download to the SSD in-repo
  `data/` (gitignored) is recommended over the HDD** — SSD is much faster for DIAMOND I/O and 55 GB
  fits in the 230 GB free. HDD remains a fine alternative if SSD space is wanted for other things.

### Hardware
- CPU: i5-12400, **12 threads** → run DIAMOND with `--threads 12`.
- GPU: RTX 3060 Ti 8 GB — **DIAMOND is CPU-only, GPU gives no speedup here**; the reference DB is
  tiny (~18k proteins) so CPU search is already fast. GPU not used by this pipeline.
- SSD /home: 230 GB free. HDD /mnt/shinrinyoku: 455 GB free (NTFS).

### gz is scientifically identical
gzip is lossless: DIAMOND (and open_text) yield byte-identical reads/bases from .gz vs .fastq.
Working from .gz changes nothing about the results — it only saves ~205 GB of disk.

## Next steps (after this)
1. Owner commits the port (see commit commands provided in chat). Owner is author; assistant not co-author.
2. Owner downloads dataset: `.venv/bin/python scripts/fetch_dataset.py` (~1.7 GB gz → ~fastq).
3. Full run: `.venv/bin/python run_pipeline.py`. Compare against `outputs_windows_baseline/`.
4. Then: research-output phase (poster/conference paper) — separate effort.

## Commit protocol
- Assistant does NOT commit. Owner (Muhammad) commits; assistant only supplies terminal commands.
- Assistant is a contributor, not commit author.
