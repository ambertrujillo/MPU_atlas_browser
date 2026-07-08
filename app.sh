#!/bin/bash
# app.sh — launch the Shiny app using the local R installation
# (native macOS setup: R packages via CRAN, Python anndata via reticulate's
# managed environment — see setup.R. No container/Singularity involved.)
set -e

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo ">>> Launching app from: $APP_DIR"
cd "$APP_DIR"

Rscript -e "shiny::runApp('app.R', launch.browser=TRUE)"