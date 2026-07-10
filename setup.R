# ---- setup.R ----
# One-time environment setup. Run this once per machine, not on every app launch.

# --- R packages ---
required_pkgs <- c("shiny", "ggplot2", "ggpubr", "reticulate", "anndata")
missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  message(">>> Installing missing R packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
} else {
  message(">>> All required R packages already installed.")
}

# --- Python side ---
library(reticulate)

# Declare the requirement up front — this is the mechanism reticulate's
# newer ephemeral-environment workflow actually expects, rather than
# relying solely on a one-off py_install() call.
reticulate::py_require("anndata")

py_ok <- tryCatch({
  py_config()
  py_module_available("anndata")
}, error = function(e) FALSE)

if (!py_ok) {
  message("!!! Python 'anndata' still not available after py_require().")
  message("!!! Manual fix: reticulate::py_install(\"anndata\")")
  message("!!! or set up a real conda environment and point RETICULATE_PYTHON at it.")
} else {
  message(">>> Python 'anndata' module found and reachable via reticulate. Setup complete.")
}

# --- Unzip datafiles ---
zip_files <- list.files("data", pattern = "\\.zip$", full.names = TRUE)

for (zip_path in zip_files) {
  expected_output <- sub("\\.zip$", "", zip_path)
  if (!file.exists(expected_output)) {
    message(">>> Unzipping ", zip_path)
    unzip(zip_path, exdir = dirname(zip_path))
  } else {
    message(">>> ", expected_output, " already exists, skipping.")
  }
}

