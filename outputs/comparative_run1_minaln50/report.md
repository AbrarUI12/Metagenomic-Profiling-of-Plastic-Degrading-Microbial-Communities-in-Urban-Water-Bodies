# Spatio-temporal profiling of plastic-degrading genetic potential in Brahmaputra river sediment

- Generated (UTC): 2026-07-20T14:52:52Z
- BioProject: PRJNA944918 (shotgun sediment metagenomes; the SRA `AMPLICON` label is a known mislabel - see Sharma et al. 2024, DOI 10.3389/fmicb.2024.1426463)

## Research question

How does the genetic potential for plastic degradation vary in space (upstream to downstream, urban versus non-urban) and across seasons in Brahmaputra river sediment, and which enzyme families and genera drive that variation?

## Study design

- Samples analysed: **18** shotgun sediment metagenomes
- Sites: **8** along the Brahmaputra
- Seasons: **3**
- Reads searched per sample: **26,000,000** (equal-depth subsampling, so hit counts are directly comparable)

| Site | Latitude | Urban | Seasons sampled |
|---|---|---|---|
| Tinsukia | 27.57 | no | 1 |
| Sadiya | 27.49 | no | 2 |
| Dibrugarh | 27.29 | no | 3 |
| Tezpur | 26.36 | no | 3 |
| Morigaon | 26.17 | no | 3 |
| Palashbari | 26.12 | no | 1 |
| Uzanbazar | 26.11 | no | 3 |
| Dhubri | 26.01 | no | 2 |

## Overall putative plastic-degrading potential

- Total filtered enzyme hits across all samples: **154,857**
- Hits by enzyme family:
  - AlkB: 78,271 (50.5%)
  - Cutinase_like: 68,490 (44.2%)
  - MHETase_like: 7,696 (5.0%)
  - PETase_like: 400 (0.3%)
- Highest potential: **Dibrugarh** (pre-monsoon (Apr 2022)), 670.46 CPM
- Lowest potential: **Tezpur** (winter / dry (Dec 2021-Jan 2022)), 87.65 CPM

## Does composition differ by site or season? (PERMANOVA)

Bray-Curtis distances on CPM profiles; pseudo-F with 999 label permutations. The balanced core (sites sampled in every season) is the primary test; the all-sample test is secondary because several sites have a single sample.

| Dataset | Factor | pseudo-F | R2 | p |
|---|---|---|---|---|
| enzyme_balanced_core | season | 50.366 | 0.918 | 0.0030 |
| enzyme_balanced_core | site | 0.032 | 0.012 | 1.0000 |
| enzyme_balanced_core | urban | 0.037 | 0.004 | 0.9070 |
| enzyme_all_samples | season | 71.406 | 0.905 | 0.0010 |
| enzyme_all_samples | site | 0.229 | 0.138 | 0.9980 |
| enzyme_all_samples | urban | 0.126 | 0.008 | 0.9280 |
| genus_balanced_core | season | 12.503 | 0.735 | 0.0010 |
| genus_balanced_core | site | 0.218 | 0.076 | 0.9960 |
| genus_balanced_core | urban | 0.258 | 0.025 | 0.9040 |
| genus_all_samples | season | 16.843 | 0.692 | 0.0010 |
| genus_all_samples | site | 0.364 | 0.203 | 1.0000 |
| genus_all_samples | urban | 0.298 | 0.018 | 0.9520 |

## Upstream -> downstream gradient

Spearman correlation against latitude. Latitude decreases downstream, so a **negative** rho means potential increases toward the downstream reaches.

| Feature | Spearman rho | p | q (BH) |
|---|---|---|---|
| total_hits_cpm | -0.041 | 0.8728 | 0.9640 |
| AlkB | -0.011 | 0.9640 | 0.9640 |
| Cutinase_like | -0.165 | 0.5140 | 0.8566 |
| MHETase_like | -0.301 | 0.2248 | 0.5619 |
| PETase_like | -0.374 | 0.1267 | 0.5619 |

## Which families and diversity metrics differ?

Kruskal-Wallis across seasons and Mann-Whitney for urban (Guwahati/Uzanbazar) versus non-urban sites, Benjamini-Hochberg corrected.

| Enzyme family | mean CPM | urban | non-urban | season p (q) | urban p (q) |
|---|---|---|---|---|---|
| Cutinase_like | 146.346 | 154.564 | 144.703 | p = 0.001 (0.002) | p = 1.000 (1.000) |
| MHETase_like | 16.444 | 18.103 | 16.113 | p = 0.002 (0.002) | p = 0.738 (1.000) |
| AlkB | 167.246 | 181.282 | 164.438 | p = 0.002 (0.002) | p = 1.000 (1.000) |
| PETase_like | 0.855 | 1.282 | 0.769 | p = 0.009 (0.009) | p = 0.155 (0.620) |

| Diversity metric | urban | non-urban | season p | urban p |
|---|---|---|---|---|
| total_hits_cpm | 355.231 | 326.023 | p = 0.002 | p = 1.000 |
| genus_richness | 477.667 | 432.667 | p = 0.002 | p = 0.498 |
| genus_shannon | 4.512 | 4.449 | p = 0.130 | p = 0.824 |
| family_shannon | 0.888 | 0.866 | p = 0.656 | p = 0.426 |

## Figures

- `outputs/comparative/figures/enzyme_composition.png`
- `outputs/comparative/figures/enzyme_heatmap.png`
- `outputs/comparative/figures/gradient_latitude.png`
- `outputs/comparative/figures/pcoa_enzyme.png`
- `outputs/comparative/figures/pcoa_genus.png`
- `outputs/comparative/figures/seasonal_boxplots.png`

## Limitations (thesis-safe)

- All results are **putative genetic potential** inferred by homology. Presence of a gene fragment is not expression, and not evidence of active plastic degradation.
- Best-hit taxonomy can misassign conserved hydrolases across taxa; genus labels are indicative, not definitive.
- The reference database is curated from public annotations and spans broad enzyme families; functional specificity requires experimental validation.
- Bulk sediment metagenomes do not demonstrate plastisphere association.
- PERMANOVA tests differences in composition; it does not establish ecological interaction, causation, or a mechanism.
- Several sites contribute a single sample, so site-level effects outside the balanced core are descriptive rather than inferential.
- Short reads (~150 bp) limit resolution; genes are detected as fragments, not assembled full-length sequences.

## Reproducibility

- Full comparative pipeline: `.venv/bin/python run_comparative.py --threads 12`
- Design table: `refs/sample_metadata.tsv`
- Tables: `outputs/comparative/tables/`, figures: `outputs/comparative/figures/`
