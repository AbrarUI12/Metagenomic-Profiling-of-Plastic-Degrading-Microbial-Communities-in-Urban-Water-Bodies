# Plastic-degrading genetic potential profiling (SRR23872596)

This project performs a thesis-grade, conservative screen for putative plastic-degrading enzymes in a polluted freshwater sediment metagenome (SRR23872596). It reports genetic potential and homology-based inference only, not active degradation.

## Quick start
1) Install Python dependencies:
```
python -m pip install -r requirements.txt
```

2) Run the full pipeline from the project root:
```
python run_pipeline.py
```

Optional quick test (subsample reads):
```
python run_pipeline.py --max-reads 200000
```

## Folder structure
- `data/` input data (optional staging)
- `refs/` reference sequences and mapping files
- `outputs/` results, tables, figures, logs, report
- `scripts/` pipeline modules
- `logs/` extra logs (reserved)
- `reports/` optional space for additional reports
- `notebooks/` optional exploration

## Outputs (key deliverables)
- `outputs/tables/genus_enzyme_matrix_raw.csv`
- `outputs/tables/genus_enzyme_matrix_cpm.csv`
- `outputs/tables/genus_clusters.csv`
- `outputs/figures/top_genera_barplot.png`
- `outputs/figures/genus_enzyme_heatmap.png`
- `outputs/figures/pca_kmeans.png`
- `outputs/figures/kmeans_elbow.png`
- `outputs/figures/kmeans_silhouette.png`
- `outputs/report.md`

## Methods overview
- Reference database: UniProtKB curated queries for PETase-like, MHETase-like, cutinase-like, AlkB, and polyesterase-like enzymes.
- Homology search: translated search (DIAMOND blastx preferred) with strict thresholds.
- Taxonomy: best-hit protein accession to genus (UniProt taxonomy).
- Clustering: k-means on log1p CPM profiles; k chosen by silhouette.

Default thresholds (configurable in `run_pipeline.py` arguments):
- e-value <= 1e-5
- percent identity >= 30
- minimum alignment length >= 50 aa

## DIAMOND installation and fallback
The pipeline will try to use DIAMOND if available. If DIAMOND is missing:
1) It attempts `conda install -y -c bioconda diamond` (if conda exists).
2) If that fails, it downloads the Windows binary from the DIAMOND GitHub release into `tools/diamond/diamond.exe`.
3) If BLAST+ is installed, it can fall back to `blastx`.
4) Otherwise, it uses a pure-Python k-mer screen (very conservative and slower, intended for small test runs).

## Reproducibility and logs
- Commands: `outputs/logs/commands.log`
- Versions: `outputs/logs/versions.txt`
- Search method log: `outputs/logs/search_method.txt`
- Parameters: `outputs/run_params.json`

## Expected runtime
- Full dataset with DIAMOND may take hours depending on CPU and I/O.
- Subsampling with `--max-reads` can complete in minutes.

## Troubleshooting
- If DIAMOND is not detected, install via conda or add to PATH.
- If the run stops at reference download, check network access.
- If memory is limited, reduce `--threads` or use `--max-reads`.

## Notes on interpretation
All results represent putative genetic potential only. Homology-based inference cannot confirm activity or plastisphere association.
