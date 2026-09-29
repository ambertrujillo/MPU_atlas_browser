#!/bin/bash

echo ">>> Starting setup for Shiny app..."

# --- Install R packages via conda ---
echo ""
echo ">>> Installing R packages via conda..."
conda install -y -c conda-forge \
  r-base \
  r-shiny \
  r-ggplot2 \
  r-ggpubr \
  r-reticulate \
  r-ggrastr \
  r-anndata

echo ">>> R packages installed via conda"

# --- Install Python packages ---
echo ""
echo ">>> Installing Python packages..."
conda install -y -c conda-forge anndata

# Install synapseclient via pip (not available in conda-forge)
pip install synapseclient

echo ">>> Python packages installed"

# --- Setup Synapse credentials ---
echo ""
echo ">>> Synapse Configuration"

SYNAPSE_CONFIG="$HOME/.synapseConfig"

if [ -f "$SYNAPSE_CONFIG" ]; then
  echo ">>> Found existing .synapseConfig at $SYNAPSE_CONFIG"
  read -p ">>> Would you like to update it? (y/n): " update_config
  
  if [ "$update_config" = "y" ]; then
    rm "$SYNAPSE_CONFIG"
    echo ">>> Removed old configuration. Creating new one..."
  else
    echo ">>> Keeping existing Synapse configuration."
  fi
fi

if [ ! -f "$SYNAPSE_CONFIG" ]; then
  echo ""
  echo ">>> Please enter your Synapse credentials:"
  echo ">>> (Get your auth token from: https://www.synapse.org -> Settings -> Personal Access Tokens)"
  
  read -p "Synapse username: " synapse_username
  read -p "Synapse auth token (Personal Access Token): " synapse_authtoken
  
  if [ -z "$synapse_username" ] || [ -z "$synapse_authtoken" ]; then
    echo "!!! Username and auth token cannot be empty."
    echo "!!! Please run setup.sh again."
    exit 1
  fi
  
  # Create .synapseConfig
  cat > "$SYNAPSE_CONFIG" << EOF
[authentication]
username = $synapse_username
authtoken = $synapse_authtoken
EOF
  
  chmod 600 "$SYNAPSE_CONFIG"
  echo ">>> Synapse configuration saved to $SYNAPSE_CONFIG"
fi

# --- Download Synapse files ---
echo ""
echo ">>> Downloading data files from Synapse..."

# Create data directory if needed
mkdir -p data

# Download files using Python
python << PYEOF
import synapseclient

syn = synapseclient.Synapse()
syn.login(silent=True)

synapse_files = [
    "syn77349469",
    "syn77615536", 
    "syn77349460"
]

for syn_id in synapse_files:
    print(f">>> Downloading {syn_id}...")
    entity = syn.get(syn_id, downloadLocation="data/")
    print(f">>> Downloaded: {entity.name}")

print(">>> All Synapse files downloaded successfully.")
PYEOF

# --- Unzip files if any ---
echo ""
echo ">>> Checking for zipped data files..."

for zip_file in data/*.zip; do
  if [ -f "$zip_file" ]; then
    basename="${zip_file%.zip}"
    if [ ! -e "$basename" ]; then
      echo ">>> Unzipping $zip_file"
      unzip -q "$zip_file" -d data/
    else
      echo ">>> $basename already exists, skipping."
    fi
  fi
done

if ! ls data/*.zip 1> /dev/null 2>&1; then
  echo ">>> No zip files found in data/ directory."
fi

echo ""
echo ">>> Setup complete! You can now run your Shiny app with: bash app.sh"
echo ">>> Synapse credentials are stored in: $SYNAPSE_CONFIG"