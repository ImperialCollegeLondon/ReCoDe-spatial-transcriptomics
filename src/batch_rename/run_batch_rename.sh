#!/bin/bash
#PBS -l walltime=0:10:0
#PBS -lselect=1:ncpus=1:mem=8gb


# Load necessary modules
module load tools/prod # Load production tools
module load Biopython/1.84-foss-2024a # Load python and bundle
source xenium_5k_venv/bin/activate # Activate virtual environment

# Batch rename script for Xenium morphology_focus files
# This script will process all output directories in your Xenium run

# Set your data directory path
DATA_DIR="/rds/general/user/sep22/projects/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/data/xenium_raw/20251103__141920__SP25164_SARA_PATTI_RUN_3"
# Get the directory where this script is located
SCRIPT_DIR="/rds/general/user/sep22/home/Projects/Xenium_5k_analysis_pipeline/src/batch_rename"

# Set your mapping file relative to this script
MAPPING_FILE="$SCRIPT_DIR/xenium_mapping.md"

echo "Xenium Batch Rename Script"
echo "=========================="
echo "Data directory: $DATA_DIR"
echo "Mapping file: $MAPPING_FILE"
echo ""

# Check if mapping file exists
if [ ! -f "$MAPPING_FILE" ]; then
    echo "Error: Mapping file '$MAPPING_FILE' not found!"
    echo "Please make sure the mapping file is in the current directory."
    exit 1
fi

echo "Proceeding with renaming..."
echo "---------------------------"
python3 "$SCRIPT_DIR/xenium_batch_rename.py" -d "$DATA_DIR" -m "$MAPPING_FILE" --rename --backup
echo ""
echo "Batch rename completed!"
echo "======================"
