# ---- setup.R ----
# One-time environment setup. Run this once per machine, not on every app launch.

# --- R packages ---
required_pkgs <- c("shiny", "ggplot2", "ggpubr", "reticulate", "anndata", "ggrastr")
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

# --- Synapse setup ---
message("\n>>> Synapse Configuration")

# Install synapseclient if needed
if (!py_module_available("synapseclient")) {
  message(">>> Installing Python synapseclient...")
  reticulate::py_require("synapseclient")
}

# Check if .synapseConfig already exists
synapse_config <- path.expand("~/.synapseConfig")

if (file.exists(synapse_config)) {
  message(">>> Found existing .synapseConfig at ", synapse_config)
  message(">>> Synapse is already configured. To reconfigure, delete ", synapse_config)
  message(">>> and run setup.R again.")
} else {
  # Credentials need to be provided via environment variables or command line args
  synapse_username <- Sys.getenv("SYNAPSE_USERNAME")
  synapse_authtoken <- Sys.getenv("SYNAPSE_AUTH_TOKEN")
  
  if (nchar(synapse_username) == 0 || nchar(synapse_authtoken) == 0) {
    message("\n>>> Synapse credentials not found.")
    message(">>> Please run setup.R with your credentials as environment variables:")
    message(">>>")
    message(">>> SYNAPSE_USERNAME='your_username' SYNAPSE_AUTH_TOKEN='your_token' Rscript setup.R")
    message(">>>")
    message(">>> Or manually create ~/.synapseConfig with:")
    message(">>> [authentication]")
    message(">>> username = your_username")
    message(">>> authtoken = your_token")
    message(">>>")
    message(">>> Get your auth token from: https://www.synapse.org -> Settings -> Personal Access Tokens")
    message("\n>>> Skipping Synapse configuration for now...")
  } else {
    # Create .synapseConfig
    config_content <- paste0(
      "[authentication]\n",
      "username = ", synapse_username, "\n",
      "authtoken = ", synapse_authtoken, "\n"
    )
    
    writeLines(config_content, synapse_config)
    
    # Set appropriate permissions (user read/write only)
    Sys.chmod(synapse_config, mode = "0600")
    
    message(">>> Synapse configuration saved to ", synapse_config)
    
    # Test Synapse login
    message("\n>>> Testing Synapse connection...")
    synapse <- tryCatch({
      syn <- import("synapseclient")
      syn$login(silent = TRUE)
      message(">>> Successfully logged in to Synapse as: ", syn$username)
      syn
    }, error = function(e) {
      message("!!! Synapse login failed: ", e$message)
      message("!!! Please verify your credentials and try running setup.R again.")
      NULL
    })
  }
}

# --- Download Synapse files ---
if (file.exists(synapse_config)) {
  message("\n>>> Downloading data files from Synapse...")
  
  # Create data directory if it doesn't exist
  if (!dir.exists("data")) {
    dir.create("data")
  }
  
  # Synapse file IDs to download
  synapse_files <- c(
    "syn77349469",  # Replace with your actual Synapse IDs
    "syn77615536",
    "syn77349460"
  )
  
  tryCatch({
    synapseclient <- import("synapseclient")
    syn <- synapseclient$Synapse()
    syn$login(silent = TRUE)
    
    for (syn_id in synapse_files) {
      message(">>> Downloading ", syn_id, "...")
      entity <- syn$get(syn_id, downloadLocation = "data/")
      message(">>> Downloaded: ", entity$name)
    }
    
    message(">>> All Synapse files downloaded successfully.")
    
  }, error = function(e) {
    message("!!! Error downloading from Synapse: ", e$message)
    message("!!! Make sure you have access to these files and your credentials are correct.")
  })
} else {
  message("\n>>> Skipping Synapse downloads (no credentials configured)")
}

# --- Unzip datafiles ---
message("\n>>> Checking for zipped data files...")
zip_files <- list.files("data", pattern = "\\.zip$", full.names = TRUE)

if (length(zip_files) == 0) {
  message(">>> No zip files found in data/ directory.")
} else {
  for (zip_path in zip_files) {
    expected_output <- sub("\\.zip$", "", zip_path)
    if (!file.exists(expected_output)) {
      message(">>> Unzipping ", zip_path)
      unzip(zip_path, exdir = dirname(zip_path))
    } else {
      message(">>> ", expected_output, " already exists, skipping.")
    }
  }
}

message("\n>>> Setup complete! You can now run your Shiny app.")
if (file.exists(synapse_config)) {
  message(">>> Synapse credentials are stored in: ", synapse_config)
}