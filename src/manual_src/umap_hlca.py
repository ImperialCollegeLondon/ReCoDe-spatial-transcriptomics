"""Filter and process the HLCA full reference dataset - Memory Optimized."""

import logging
import sys
from pathlib import Path

import scanpy as sc

log_file = "hlca_analysis.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_file),
    ],
)

logging.info("Starting HLCA plotting script...")

# Set directories and file paths
output_dir = Path("/rds/general/user/sep22/home/Projects/_Public_datasets/HLCA/data/")
adata_path = output_dir / "hlca_full_unprocessed.h5ad"
filtered_path = output_dir / "hlca_full_filtered.h5ad"
processed_path = output_dir / "hlca_full_processed.h5ad"

# Set figure directory for this module
sc.settings.figdir = output_dir / "figs" / "hlca_full_filt_process"
sc.settings.figdir.mkdir(parents=True, exist_ok=True)

# Check if adata file exists
if not adata_path.exists():
    raise FileNotFoundError(f"{adata_path} not found.")

# Use backed mode to avoid loading entire dataset into memory
adata = sc.read_h5ad(processed_path)

# Plot UMAP
obs_lists = ["disease", "tissue_level_2", "cell_type"]
for obs in obs_lists:
    sc.pl.umap(
        adata,
        color=obs,
        title=f"HLCA full reference data UMAP colored by {obs}",
        save=f"_hlca_full_ref_{obs}.pdf",
    )
logging.info("UMAP plots saved.")
