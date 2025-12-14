#!/bin/bash
#PBS -l walltime=8:0:0
#PBS -lselect=1:ncpus=32:mem=512gb:ngpus=1:gpu_type=A100
#PBS -q v1_gpu72

# Load production tools
module load tools/prod

# Load python and bundle
module load Biopython/1.84-foss-2024a

# Load CUDA if needed
module load CUDA/12.1.1

# Change to directory
cd /rds/general/user/sep22/home/Projects/Xenium_5k_analysis_pipeline

# Activate virtual environment
source xenium_5k_venv/bin/activate

# Set GPU environment variables
export CUDA_VISIBLE_DEVICES=0

# Run
python -m recode_st config_files/config_HPC_runs_int.toml
