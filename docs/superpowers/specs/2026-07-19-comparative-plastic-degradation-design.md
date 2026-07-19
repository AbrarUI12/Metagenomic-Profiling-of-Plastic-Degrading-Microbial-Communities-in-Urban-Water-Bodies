# Design: Spatio-temporal comparative screen of plastic-degrading genetic potential across the Brahmaputra

- Date: 2026-07-19
- Status: Draft for review
- Scope: One implementation plan (adds a comparative analysis layer; does not modify the existing single-sample pipeline)

## 1. Motivation and novelty

The existing project screens **one** sediment metagenome (SRR23872596, Tezpur) for
homologs of plastic-degrading enzymes and reports genetic potential for that single
sample. We now hold the **entire BioProject PRJNA944918 — 18 runs** covering **8 sites
along ~600 km of the Brahmaputra, sampled in 3 seasons**.

The published source study (Sharma et al. 2024, Front. Microbiol., DOI
10.3389/fmicb.2024.1426463) describes only general microbial community composition,
functional annotation, and antibiotic-resistance genes for the 6 post-monsoon samples. It
does **not** address plastic-degrading enzymes, nor the spatial/seasonal structure of any
functional trait. Our contribution:

> **How does the genetic potential for plastic degradation vary in space (upstream →
> downstream, urban vs non-urban) and season (post-monsoon / winter / pre-monsoon) in
> Brahmaputra river sediment, and which enzyme families and genera drive that variation?**

This turns a single-sample "genetic potential" screen into a statistically testable
spatio-temporal comparison — the natural, publishable upgrade.

## 2. Data and experimental design

18 runs, all "sediment metagenome", `library_strategy=AMPLICON` in SRA but confirmed
**shotgun** (see `handout.md`; cite the mislabel resolution). Structure derived from ENA
metadata (`sample_title`, `lat`, `lon`, `collection_date`):

**8 sites (upstream → downstream by latitude / river position):**

| Site | Code | lat, lon | Note |
|---|---|---|---|
| Tinsukia | TSB | 27.57, 95.32 | most upstream (near source) |
| Sadiya | SAB/SBR | 27.49, 95.40 | upstream |
| Dibrugarh | DGB | 27.29, 94.54 | upstream town |
| Tezpur | TZB | 26.36, 92.47 | mid (original sample) |
| Morigaon | MRB | 26.17, 92.06 | mid |
| Uzanbazar | UZB | 26.11, 91.45 | **Guwahati — largest city on the river (~1.1 M); the "urban water body"** |
| Palashbari | PAB | 26.12, 91.54 | mid-downstream |
| Dhubri | DHB | 26.01, 89.59 | most downstream (near Bangladesh border) |

**3 seasonal campaigns (by collection date):**
- **Post-monsoon** — Sep 2021 (6 samples; the original published BRS 1–6)
- **Winter / dry** — Dec 2021 – Jan 2022 (8 samples)
- **Pre-monsoon** — Apr 2022 (4 samples)

**Balanced core:** Uzanbazar, Morigaon, Tezpur, Dibrugarh each have all 3 seasons →
**4 sites × 3 seasons = 12 samples** fully crossed. The remaining 6 samples add spatial
coverage (Dhubri ×2, Sadiya ×2, Palashbari ×1, Tinsukia ×1).

Season labels will be assigned from `collection_date` (ground truth), not from the SRA
sample-code suffixes, which are ambiguous.

## 3. Key methodological decisions

- **Depth normalization (decided): subsample every sample to a common depth of
  ~13 M read pairs** (≈ the smallest sample, Tezpur post-monsoon at 13.5 M) before the
  homology search. Sequencing depth ranges 13 M–70 M pairs; without equalization, deeper
  samples yield more enzyme hits purely as an artifact. Equal-depth subsampling makes
  hit counts directly comparable and also bounds compute. CPM is still reported as a
  secondary normalization. Subsampling is deterministic (fixed seed) for reproducibility.
- **No new heavy dependencies.** `skbio`/`statsmodels` are not installed. PCoA (classical
  MDS) and PERMANOVA are implemented in ~30 lines each on top of `numpy`/`scipy`
  (`scipy.spatial.distance` for Bray–Curtis, `numpy.linalg.eigh` for PCoA). This keeps the
  environment reproducible and portable, consistent with the project's existing choices.
- **Existing single-sample pipeline is untouched.** `run_pipeline.py` and `scripts/00..99`
  keep working exactly as before. The comparative work is a parallel layer.
- **Reuse the existing reference DB and search logic.** The same UniProt enzyme FASTA
  (`refs/plastic_enzymes.fasta`) and DIAMOND blastx thresholds (e-value ≤ 1e-5, pident
  ≥ 30, min aln 50 aa) are applied per sample, so per-sample results stay comparable to the
  original single-sample report.

## 4. Analysis architecture

New orchestrator `run_comparative.py` runs a new numbered series `10 → 14`. Data flows
one direction; each stage has a single responsibility and a file interface.

```
refs/sample_metadata.tsv        (site, code, lat, lon, season, date, urban flag, run)
        │
        ▼
10_batch_search.py   per sample: subsample→DIAMOND→taxonomy   (RESUMABLE, skips done)
        │   outputs/comparative/per_sample/<RUN>/enzyme_hits_with_genus.tsv
        ▼
11_build_matrices.py   combine 18 samples
        │   outputs/comparative/tables/sample_enzyme_matrix_{raw,cpm}.csv
        │   outputs/comparative/tables/sample_genus_matrix_{raw,cpm}.csv
        │   outputs/comparative/tables/sample_summary.tsv   (per-sample totals + metadata)
        ▼
12_comparative_stats.py   PCoA, PERMANOVA, diversity, differential abundance, gradient
        │   outputs/comparative/tables/*.csv  (stats results)
        ▼
13_comparative_plots.py   ordination / gradient / seasonal / heatmap figures
        │   outputs/comparative/figures/*.png
        ▼
14_comparative_report.py   outputs/comparative/report.md  (results-section style)
```

### 4.1 `refs/sample_metadata.tsv` (committed)
18 rows built once from the ENA pull. Columns: `run_accession, site, site_code, lat, lon,
season, collection_date, urban (0/1), read_pairs_total`. `urban=1` only for Uzanbazar
(Guwahati). Season from date. Committed so the analysis is deterministic and the ENA call
is not repeated on every run. A small `scripts/make_sample_metadata.py` generates it (and
can be re-run to refresh), but the committed TSV is the source of truth for the analysis.

### 4.2 `scripts/10_batch_search.py`
For each run in the metadata table:
1. Subsample R1/R2 to N read pairs (default 13 M, `--depth`) deterministically. Prefer
   `seqtk sample -s<seed>` if available; otherwise a pure-Python paired reservoir/streaming
   sampler over the gzipped FASTQ (reads `.fastq.gz` directly via existing `open_text`).
2. DIAMOND blastx of the subsample vs `refs/plastic_enzymes.fasta` with the standard
   thresholds (reuse `03_diamond_search.py` internals or invoke it with per-sample paths).
3. Best-hit taxonomy → genus (reuse `04_taxonomy_map.py` logic).
4. Write `outputs/comparative/per_sample/<RUN>/enzyme_hits_with_genus.tsv` and a tiny
   `<RUN>/sample_meta.json` (depth used, hit count, seed).

**Resumable:** a sample whose output exists and is non-empty is skipped. This mirrors
`fetch_dataset.py`'s crash-safety so an interrupted batch resumes cleanly. Writes to a
`.part` then renames.

**Important:** must NOT write to the single-sample hardcoded paths under
`outputs/enzyme_hits/` (the known clobber gotcha in `handout.md`). All comparative output
lives under `outputs/comparative/`.

### 4.3 `scripts/11_build_matrices.py`
Reads the 18 per-sample hit tables. Produces:
- `sample_enzyme_matrix_raw.csv` / `_cpm.csv` — rows = samples, cols = enzyme families
  (AlkB, Cutinase_like, PETase_like, MHETase_like, Polyesterase). CPM = hits / subsample
  depth × 1e6.
- `sample_genus_matrix_raw.csv` / `_cpm.csv` — rows = samples, cols = genera (long tail;
  keep all, downstream code selects top-N for plots).
- `sample_summary.tsv` — per sample: total hits, hits per family, Shannon diversity of
  genera and of enzyme families, plus the metadata columns (site, season, lat, urban).

### 4.4 `scripts/12_comparative_stats.py`
- **Ordination:** Bray–Curtis distance on `sample_enzyme_matrix_cpm` and on
  `sample_genus_matrix_cpm`; classical PCoA (numpy eigh). Save coordinates +
  variance-explained per axis.
- **PERMANOVA** (pseudo-F, 999 permutations, fixed seed) on the Bray–Curtis distances,
  testing `~ season` and `~ site` (and `~ urban`). Report pseudo-F, R², p. Run the
  primary two-factor test on the **balanced 4×3 core** (cleanest design); report the
  full-18 test as secondary.
- **Alpha diversity** (Shannon, richness) per sample; compare across season and across
  urban/non-urban with a nonparametric test (Kruskal–Wallis via `scipy.stats`).
- **Differential abundance:** per enzyme family and per top genus, compare urban vs
  non-urban and season groups (Kruskal–Wallis / Mann–Whitney, Benjamini–Hochberg FDR).
- **Spatial gradient:** Spearman correlation of per-sample total plastic-potential (CPM)
  and per-family CPM against latitude (upstream→downstream proxy).
- All results written as CSVs; nothing hardcoded assumes a specific significant result.

### 4.5 `scripts/13_comparative_plots.py`
- PCoA scatter, points colored by **site**, shaped by **season** (functional and taxonomic
  versions).
- Latitude (upstream→downstream) vs total potential, with per-season trend.
- Seasonal boxplots of total potential and of each enzyme family (balanced core).
- Sample × enzyme-family heatmap (CPM, clustered), annotated by site/season.
- Stacked bar of enzyme-family composition per sample, grouped by site then season.

### 4.6 `scripts/14_comparative_report.py`
Markdown report in results-section style: data/design table, ordination + PERMANOVA
outcome, gradient result, seasonal result, differential families/genera, figure list, and
the carried-over thesis-safe limitations. Numbers are read from the stats CSVs, never
hardcoded.

## 5. Deliverables

- `outputs/comparative/` — matrices, stats tables, figures, `report.md`.
- `refs/sample_metadata.tsv` — the committed design table.
- The figure + stats set is structured so it can drop directly into a poster or a paper's
  Results section later. The output *format* (poster vs paper vs thesis chapter) is a
  separate downstream choice and is out of scope for this plan.

## 6. Thesis-safe framing (carried over, unchanged)

All results are **putative genetic potential** by homology-based inference — presence of a
gene is not expression or activity; best-hit taxonomy can misassign conserved hydrolases;
bulk sediment does not prove plastisphere association. Spatial/seasonal *differences in
potential* inherit every one of these limitations and must be stated as such. PERMANOVA
tests differences in composition, not ecological interaction or function.

## 7. Out of scope (YAGNI)

- No metatranscriptomics/metaproteomics, no assembly/MAGs, no wet-lab (these remain
  "future directions" in `explanation.md`).
- No change to the single-sample pipeline or its committed outputs.
- No new enzyme families or reference DB changes — same DB as the baseline, so results are
  comparable.
- No PDF/poster rendering in this plan — deliver analysis + markdown; formatting is a later
  step.

## 8. Open questions / risks

- **Subsample sampler:** `seqtk` is not installed. The pure-Python paired sampler must be
  correct (keep R1/R2 in sync) and memory-safe over 13 M pairs. Plan should unit-test it on
  a tiny synthetic pair. (Adding `seqtk` via a static binary is an alternative if the Python
  sampler is too slow.)
- **Runtime:** 18 × DIAMOND blastx on 13 M-pair subsamples vs an ~18k-protein DB on 12
  threads. Expected minutes–tens-of-minutes per sample; full batch a few hours. Resumability
  makes this safe to run in stages.
- **Few replicates per site** (some sites have 1 sample). Mitigated by anchoring the primary
  statistical test on the balanced 4×3 core and treating single-sample sites as spatial
  context, not as their own test group.
