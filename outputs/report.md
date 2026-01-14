# Plastic-degrading genetic potential profiling report
- Generated (UTC): 2026-01-14T11:36:28.379638Z

## Research topic
Using shotgun metagenomic reads from a polluted freshwater environment, screen for homologs of known plastic-degrading enzymes (e.g., PETase/cutinase-like hydrolases, MHETase-like, alkane monooxygenase AlkB, relevant esterases/lipases where appropriate). Assign each enzyme hit to a likely microbial genus by best-hit taxonomy (homology), quantify enzyme-hit abundance per genus, then perform k-means clustering of genera based on enzyme-hit profiles to identify 'high / moderate / low plastic-degradation potential' groups. Report results with conservative, thesis-safe language emphasizing genetic potential and limitations.

## Data description
- Accession: SRR23872596
- Total reads (R1+R2): 26,990,268
- Read1 length (min/mean/max): 159 / 159.00 / 159
- Read2 length (min/mean/max): 159 / 159.00 / 159
- Key metadata fields:
  - BioProject: PRJNA944918
  - BioSample: SAMN33764824
  - LibrarySelection: PCR
  - LibrarySource: METAGENOMIC
  - LibraryStrategy: AMPLICON
  - LoadDate: 2023-03-17 12:19:24
  - Platform: ILLUMINA
  - ReleaseDate: 2024-05-01 07:18:45
  - Run: SRR23872596
  - Sample: SRS17053977
  - SampleName: TezpurPTM_Brahmaputra
  - ScientificName: sediment metagenome
  - avgLength: 318
  - bases: 4291452612
  - spots: 13495134

## Methods (homology-based inference)
- Reference database: UniProtKB stream queries for PETase-like, MHETase-like, cutinase-like, AlkB, and polyesterase-like enzymes (see `refs/README_refs.md`).
- Enzyme detection: translated search of reads vs protein database; top hit per read retained.
- Thresholds: e-value <= 1e-05, percent identity >= 30.0, min alignment length >= 50.
- Subsampling: analysis limited to first 2000000 reads per FASTQ for this run.
- Search method log: `outputs/logs/search_method.txt`
- Taxonomy assignment: best-hit mapping from UniProt protein accessions to genus (homology-based).
- Clustering: k-means on log1p CPM enzyme-family profiles; k selected by silhouette.

## Results (putative genetic potential)
- Total putative enzyme hits: 2,060 across 277 genera.
- Top genera by total putative plastic-enzyme hits:
  - Aspergillus: 266
  - Alcanivorax: 113
  - Fusarium: 95
  - Gordonia: 65
  - Streptomyces: 63
  - Pyrenophora: 56
  - Sphingopyxis: 56
  - Mycobacterium: 55
  - Acaromyces: 51
  - Paraperlucidibaca: 47
- Top enzyme families by hit count:
  - AlkB: 1,134
  - Cutinase_like: 839
  - MHETase_like: 81
  - PETase_like: 6

Cluster summary (putative plastic-degradation potential groups):
- Cluster 0: 273 genera
- Cluster 1: 4 genera

## Network-Based Modularity Analysis of Plastic-Degrading Potential
- Rationale: graph modularity can reveal putative functional modules (guilds) based on co-occurrence of plastic-degrading genetic potential, complementing distance-based clustering.
- Network construction: genus-genus similarity computed from CPM enzyme-family profiles using cosine with threshold 0.6; edge weights retained and self-loops removed.
- Community detection: louvain (graph modularity), seed 42.
- Modularity (Q): 0.255; modules: 4; nodes: 277; edges: 20613.
- Major modules (by total hits):
  - Module 0: 179 genera, 1,137 total hits, top families: AlkB;Cutinase_like;MHETase_like
  - Module 1: 82 genera, 853 total hits, top families: Cutinase_like;AlkB;MHETase_like
  - Module 2: 12 genera, 62 total hits, top families: MHETase_like;AlkB
  - Module 3: 4 genera, 8 total hits, top families: PETase_like;AlkB;MHETase_like
- Comparison with k-means: overlap table saved to `outputs/tables/kmeans_vs_modules.csv`; alignment is partial and reflects complementary views of functional organization.
- Adjusted Rand Index (modules vs k-means): 0.058
- Interpretation: modules represent putative functional communities structured by shared plastic-degradation genetic potential (homology-based inference), not confirmed interactions or activity.

## Figures
- `outputs/figures/top_genera_barplot.png`
- `outputs/figures/genus_enzyme_heatmap.png`
- `outputs/figures/pca_kmeans.png`
- `outputs/figures/kmeans_elbow.png`
- `outputs/figures/kmeans_silhouette.png`
- `outputs/figures/network/genus_network_modules.png`
- `outputs/figures/network/genus_enzyme_heatmap_by_module.png`
- `outputs/figures/network/kmeans_vs_modules_heatmap.png`
- `outputs/figures/network/module_barplots.png`

## Limitations (thesis-safe)
- These are putative, candidate plastic-degrading genes identified by homology-based inference; metagenomic presence does not imply expression or active degradation.
- Best-hit taxonomy can misassign genes across taxa, especially for conserved hydrolases.
- The database is curated from public annotations and may include broad enzyme families; functional specificity requires experimental validation.
- Results describe genetic potential in bulk sediment metagenome and do not prove plastisphere association.
- Correlation-based networks reflect co-occurrence of genetic potential, not ecological interaction.
- Network structure and modularity depend on similarity metric and threshold choices.

## Reproducibility
- Run end-to-end pipeline: `python run_pipeline.py`
- Outputs are saved under `outputs/` with logs in `outputs/logs/`.
