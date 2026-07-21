# Spatio-temporal profiling of plastic-degrading genetic potential in Brahmaputra river sediment

- Generated (UTC): 2026-07-20T16:49:16Z
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

- Total filtered enzyme hits across all samples: **441,298**
- Hits by enzyme family:
  - Cutinase_like: 222,555 (50.4%)
  - AlkB: 197,440 (44.7%)
  - MHETase_like: 19,834 (4.5%)
  - PETase_like: 1,466 (0.3%)
  - Polyesterase_like: 3 (0.0%)
- Highest potential: **Dibrugarh** (pre-monsoon (Apr 2022)), 1156.81 CPM
- Lowest potential: **Tezpur** (winter / dry (Dec 2021-Jan 2022)), 684.88 CPM

## Does composition differ by site or season? (PERMANOVA)

Bray-Curtis distances on CPM profiles; pseudo-F with 999 label permutations. The balanced core (sites sampled in every season) is the primary test; the all-sample test is secondary because several sites have a single sample.

| Dataset | Factor | pseudo-F | R2 | p |
|---|---|---|---|---|
| enzyme_balanced_core | season | 3.058 | 0.405 | 0.0480 |
| enzyme_balanced_core | site | 0.361 | 0.119 | 0.8930 |
| enzyme_balanced_core | urban | 0.054 | 0.005 | 0.9230 |
| enzyme_all_samples | season | 0.980 | 0.116 | 0.3920 |
| enzyme_all_samples | site | 0.425 | 0.229 | 0.9330 |
| enzyme_all_samples | urban | 0.000 | 0.000 | 1.0000 |
| genus_balanced_core | season | 2.349 | 0.343 | 0.0040 |
| genus_balanced_core | site | 0.686 | 0.205 | 0.8740 |
| genus_balanced_core | urban | 0.627 | 0.059 | 0.8040 |
| genus_all_samples | season | 1.768 | 0.191 | 0.0280 |
| genus_all_samples | site | 0.705 | 0.330 | 0.9280 |
| genus_all_samples | urban | 0.570 | 0.034 | 0.9070 |

## Upstream -> downstream gradient

Spearman correlation against latitude. Latitude decreases downstream, so a **negative** rho means potential increases toward the downstream reaches.

| Feature | Spearman rho | p | q (BH) |
|---|---|---|---|
| total_hits_cpm | 0.104 | 0.6808 | 0.9020 |
| AlkB | 0.205 | 0.4140 | 0.8280 |
| Cutinase_like | 0.080 | 0.7517 | 0.9020 |
| MHETase_like | -0.406 | 0.0943 | 0.5660 |
| PETase_like | -0.274 | 0.2708 | 0.8123 |
| Polyesterase_like | 0.015 | 0.9545 | 0.9545 |

## Which families and diversity metrics differ?

Kruskal-Wallis across seasons and Mann-Whitney for urban (Guwahati/Uzanbazar) versus non-urban sites, Benjamini-Hochberg corrected.

| Enzyme family | mean CPM | urban | non-urban | season p (q) | urban p (q) |
|---|---|---|---|---|---|
| Cutinase_like | 475.545 | 476.244 | 475.405 | p = 0.079 (0.315) | p = 0.738 (0.912) |
| PETase_like | 3.132 | 3.744 | 3.010 | p = 0.126 (0.315) | p = 0.286 (0.912) |
| MHETase_like | 42.380 | 44.910 | 41.874 | p = 0.282 (0.470) | p = 0.574 (0.912) |
| Polyesterase_like | 0.006 | 0.013 | 0.005 | p = 0.427 (0.534) | p = 0.464 (0.912) |
| AlkB | 421.880 | 413.282 | 423.600 | p = 0.622 (0.622) | p = 0.912 (0.912) |

| Diversity metric | urban | non-urban | season p | urban p |
|---|---|---|---|---|
| total_hits_cpm | 938.192 | 943.895 | p = 0.725 | p = 1.000 |
| genus_richness | 673.333 | 643.733 | p = 0.373 | p = 0.574 |
| genus_shannon | 4.545 | 4.494 | p = 0.120 | p = 0.912 |
| family_shannon | 0.868 | 0.852 | p = 0.199 | p = 0.426 |

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
