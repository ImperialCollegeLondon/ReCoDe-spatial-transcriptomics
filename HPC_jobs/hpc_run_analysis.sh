#!/bin/bash
#PBS -l walltime=1:0:0
#PBS -lselect=1:ncpus=8:mem=256gb

# Load production tools
module load tools/prod

# Load python and bundle
module load Biopython/1.84-foss-2024a

# Change to directory
cd /rds/general/user/sep22/home/Projects/Xenium_5k_analysis_pipeline

# Activate virtual environment
source xenium_5k_venv/bin/activate

# Run
python -m recode_st config_files/config_HPC_drug2cell.toml
