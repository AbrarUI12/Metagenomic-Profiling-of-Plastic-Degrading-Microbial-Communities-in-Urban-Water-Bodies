# Plastic-degrading genetic potential profiling report
- Generated (UTC): 2026-01-14T08:26:10.377546Z

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
- Subsampling: analysis limited to first 200000 reads per FASTQ for this run.
- Search method log: `outputs/logs/search_method.txt`
- Taxonomy assignment: best-hit mapping from UniProt protein accessions to genus (homology-based).
- Clustering: k-means on log1p CPM enzyme-family profiles; k selected by silhouette.

## Results (putative genetic potential)
- Total putative enzyme hits: 212 across 78 genera.
- Top genera by total putative plastic-enzyme hits:
  - Alcanivorax: 18
  - Aspergillus: 17
  - Fusarium: 11
  - Pyrenophora: 10
  - Streptomyces: 10
  - Gordonia: 8
  - Alicyclobacillus: 6
  - Mycobacterium: 6
  - marine: 5
  - Paraperlucidibaca: 5
- Top enzyme families by hit count:
  - AlkB: 115
  - Cutinase_like: 87
  - MHETase_like: 10

Cluster summary (putative plastic-degradation potential groups):
- Cluster 0: 72 genera
- Cluster 1: 3 genera
- Cluster 2: 3 genera

## Figures
- `outputs/figures/top_genera_barplot.png`
- `outputs/figures/genus_enzyme_heatmap.png`
- `outputs/figures/pca_kmeans.png`
- `outputs/figures/kmeans_elbow.png`
- `outputs/figures/kmeans_silhouette.png`

## Limitations (thesis-safe)
- These are putative, candidate plastic-degrading genes identified by homology-based inference; metagenomic presence does not imply expression or active degradation.
- Best-hit taxonomy can misassign genes across taxa, especially for conserved hydrolases.
- The database is curated from public annotations and may include broad enzyme families; functional specificity requires experimental validation.
- Results describe genetic potential in bulk sediment metagenome and do not prove plastisphere association.

## Reproducibility
- Run end-to-end pipeline: `python run_pipeline.py`
- Outputs are saved under `outputs/` with logs in `outputs/logs/`.
