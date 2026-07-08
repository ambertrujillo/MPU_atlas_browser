# CARD-MPU Atlas Browser
 
An interactive Shiny app for exploring single-nucleus RNA-seq data from the CARD-MPU atlas', spanning Alzheimer's disease (AD) and frontotemporal dementia (FTD) across the Trujillod, Marsan, and Mathys cohorts.
 
The app provides two main views:
 
- **Pseudobulk expression browser** — violin plots of per-sample gene expression across disease groups (with optional pairwise statistical comparisons), and cell-type barplots of expression by disease group.
- **Feature plot browser** — single-cell UMAP plots colored by gene expression, drawn from the full thalamus atlas or a neuron-subtype-only object.
---
 
## App Setup
 
### Prerequisites
- R (version 4.4 or later recommended)
- Internet access on first run (for package installation)
### 1. Install dependencies
Run this once, the first time you set up the app on a new machine:
```bash
Rscript setup.R
```
This installs the required R packages (`shiny`, `ggplot2`, `ggpubr`, `reticulate`, `anndata`) and configures the Python environment used for AnnData/single-cell feature plots. You only need to run this once — subsequent launches will skip anything already installed.
 
### 2. Launch the app
```bash
bash app.sh
```
This starts the Shiny app and opens it in your default browser.
 
---
 
## Data requirements
 
The app expects the following files to be present locally under `data/`:
 
| File | Used by |
|---|---|
| `Thalamus_SampleID.csv` | Sample-level pseudobulk violin plot |
| `Thalamus_generalDisease_celltype.csv` | Cell-type barplot |
| `Thalamus_generalDisease_celltype_splitFTD` | Cell-type barplot, split by FTD type |
| `combined_reclustered.h5ad` | Feature plot — full atlas |
| `Neu_adata_annotation_04142026.h5ad` | Feature plot — neuron subtypes only |
 
These files are not included in this repository due to size. Contact [your name / group] for access, and place them in `data/` before launching the app.
 
---
 
## Troubleshooting
 
**`setup.R` reports Python `anndata` not found** — re-run `Rscript setup.R`; it will attempt to install it automatically into the environment `reticulate` is bound to. If the issue persists, see the comments in `setup.R` for manual fix options (e.g., `reticulate::py_install("anndata")` or setting up a named conda environment).
 
**"File not found" errors when loading data** — confirm the files listed above are present under `data/` with the exact filenames shown.
 
**Feature plot fails after a successful AnnData load** — confirm the embedding key (`X_umap_mnn`) and expression layer (`cpm`) referenced in `app.R` actually exist in the loaded object; the app will report which keys are available if there's a mismatch.
 
---
 
## Project structure
```
.
├── app.R          # Main Shiny application
├── setup.R        # One-time dependency installation
├── app.sh         # App launcher
├── data/          # Input CSVs and h5ad files (not tracked in git)
└── www/           # Static assets (logo, reference images)
```
 