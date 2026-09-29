![screenshot](www/card_mpu_logo.svg)
 
An interactive Shiny app for exploring single-nucleus RNA-seq data from the CARD-MPU atlas', spanning Alzheimer's disease (AD) and frontotemporal dementia (FTD) across the Trujillo, Marsan, and Mathys cohorts.
 
The app provides two main views:
 
- **Pseudobulk expression browser** — violin plots of per-sample gene expression across disease groups (with optional pairwise statistical comparisons), and cell-type barplots of expression by disease group.
- **Feature plot browser** — single-cell UMAP plots colored by gene expression, drawn from the full thalamus atlas or a neuron-subtype-only object.
---
 
## App Setup
 
### Prerequisites
- R (version 4.4 or later recommended)
```bash
conda install -c conda-forge r-base
```
- Internet access on first run (for package installation)

### 1. Install dependencies
Run this once, interactively, the first time you set up the app on a new machine:
```bash
bash setup.sh
```
You will be prompted to enter your Synapse username: <your_username> and Synapse auth token (Personal Access Token): <your_access_token>.

This installs the required R packages (`shiny`, `ggplot2`, `ggpubr`, `reticulate`, `anndata`, `ggrastr`) and configures the Python environment used for AnnData/single-cell feature plots. You only need to run this once — subsequent launches will skip anything already installed.
 
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
| `combined_reclustered_trimmed.h5ad` | Feature plot — full atlas |
| `Neu_thalamus_adata.h5ad` | Feature plot — neuron subtypes only |
 
These files are not included in this repository due to size. Contact [your name / group] for access, and place them in `data/` before launching the app.
 
---
 
## Troubleshooting

**Python environment issues (anndata not found)** — The app uses `reticulate` to access Python's anndata library. If you see "ModuleNotFoundError: No module named 'anndata'", reticulate may be using a different Python than where conda installed packages. Fix by adding this line to the very top of `app.R` (before any library calls):
```r
Sys.setenv(RETICULATE_PYTHON = "/opt/miniconda3/bin/python")
```
Then restart your R session. Verify with `reticulate::py_config()` and `reticulate::py_module_available("anndata")`.

**Missing R or Python packages after setup** — Re-run `./setup.sh`. It will install any missing conda packages (R and Python) and configure Synapse credentials if needed.

**"File not found" errors when loading data** — Confirm that `setup.sh` successfully downloaded files to `data/` and that the filenames in `app.R` match exactly. Check:
- `data/combined_reclustered_trimmed.h5ad`
- `data/Neu_thalamus_adata.h5ad`  
- `data/Thalamus_*.csv` files

**Feature plot fails after successful AnnData load** — Verify that the embedding key (`X_umap_mnn`) and expression layer (`cpm`) exist in your h5ad file. The app will report available keys if there's a mismatch. Check your AnnData object structure matches what the app expects.

**Synapse download fails during setup** — Verify your Personal Access Token has "Download" permissions enabled (not just "View"). Create a new token at https://www.synapse.org → Settings → Personal Access Tokens, making sure to check the "Download" scope.

**Cairo/XQuartz errors** — The app uses standard PDF output and doesn't require Cairo. If you see Cairo-related errors, they should not prevent the app from running. Graphics will still render correctly.
 
---
 
## Project structure
```
.
├── app.R          # Main Shiny application
├── setup.sh        # One-time dependency installation
├── app.sh         # App launcher
├── data/          # Input CSVs and h5ad files (not tracked in git)
└── www/           # Static assets (logo, reference images)
```
 