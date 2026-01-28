"""Filter and process the HLCA full reference dataset - Memory Optimized."""

import logging
import sys
from pathlib import Path

import numpy as np
import scanpy as sc

# Global variables
REF_DATASET = "hlca_full_ref"

# Set directories and file paths
output_dir = Path("/rds/general/user/sep22/home/Projects/_Public_datasets/HLCA/data/")
adata_path = output_dir / f"unprocessed_{REF_DATASET}.h5ad"
filtered_path = output_dir / f"filtered_{REF_DATASET}.h5ad"
processed_path = output_dir / f"processed_{REF_DATASET}.h5ad"


# Define a functions
def verify_counts_layer(adata, context=""):
    """Verify that counts layer exists and contains valid count data."""
    if "counts" not in adata.layers:
        logging.error(f"CRITICAL: Counts layer missing {context}!")
        raise ValueError(f"Counts layer was lost {context}")

    logging.info(f"✓ Counts layer verified: shape {adata.layers['counts'].shape}")

    # Verify it contains actual count data
    counts_sample = adata.layers["counts"][:100, :100]
    counts_arr = counts_sample.A if hasattr(counts_sample, "A") else counts_sample
    is_integer = np.allclose(counts_arr, np.round(counts_arr))
    has_counts = np.median(counts_arr.sum(axis=1)) > 100

    if is_integer and has_counts:
        logging.info("✓ Counts layer contains valid count data")
        return True
    else:
        logging.warning("⚠ Counts layer may not contain valid count data")
        return False


# Define logging configuration
log_file = output_dir / f"{REF_DATASET}_analysis.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_file),
    ],
)

logging.info("Starting reference plotting script...")
logging.info(f"Reference dataset: {REF_DATASET}")


# Set figure directory for this module
sc.settings.figdir = output_dir / "figs" / REF_DATASET
sc.settings.figdir.mkdir(parents=True, exist_ok=True)

# Check if adata file exists
if not adata_path.exists():
    raise FileNotFoundError(f"{adata_path} not found.")

logging.info(f"Loading adata file from {adata_path}...")
# Use backed mode to avoid loading entire dataset into memory
adata = sc.read_h5ad(adata_path, backed="r")
logging.info(f"Original shape: {adata.shape}")

logging.info("Filter for disease and tissue_level_2...")
for key in ["disease", "tissue_level_2", "cell_type"]:
    if key not in adata.obs:
        raise KeyError(f"{key} not found in adata.obs columns.")

# Create filter mask
keep_diseases = [
    "normal",
    "pulmonary fibrosis",
    "interstitial lung disease",
    "chronic obstructive pulmonary disease",
]
remove_tissues = ["nose", "inferior turbinate", "trachea"]

disease_mask = adata.obs["disease"].isin(keep_diseases)
tissue_mask = ~adata.obs["tissue_level_2"].isin(remove_tissues)
combined_mask = disease_mask & tissue_mask

logging.info(f"Cells passing filters: {combined_mask.sum()} / {len(combined_mask)}")

# Load only filtered subset into memory
logging.info("Loading filtered subset into memory...")
adata = adata[combined_mask, :].to_memory()
logging.info(f"Filtered shape: {adata.shape}")

logging.info(f"Remaining diseases: {adata.obs['disease'].unique()}")
logging.info(f"Remaining tissue types: {adata.obs['tissue_level_2'].unique()}")

# Confirm filtering
logging.info(f"Disease value counts:\n{adata.obs['disease'].value_counts()}")
logging.info(f"Tissue value counts:\n{adata.obs['tissue_level_2'].value_counts()}")

# Save filtered dataset
logging.info(f"Saving filtered adata to {filtered_path}...")
adata.write_h5ad(filtered_path)

logging.info("Process filtered dataset...")

# Store raw counts BEFORE any filtering
logging.info("Checking for counts layer before filtering...")
if "counts" not in adata.layers:
    logging.info("Counts layer not found - attempting to create from adata.raw...")
    if adata.raw is not None:
        logging.info("✓ adata.raw found - creating counts layer from raw data")
        adata.layers["counts"] = adata.raw.X.copy()
        logging.info("✓ Successfully created counts layer from adata.raw")
    else:
        logging.error("✗ No adata.raw found - cannot create counts layer")
        logging.error("Please ensure your data has either:")
        logging.error("  1. An existing 'counts' layer with raw count data, or")
        logging.error("  2. An adata.raw attribute with the original count matrix")
        raise ValueError(
            "Cannot create counts layer: no 'counts' layer found and "
            "no adata.raw available"
        )
else:
    logging.info("✓ Counts layer found - using existing counts layer")

# Verify the counts layer before filtering
verify_counts_layer(adata, "before filtering")

# Basic filtering
logging.info("Filtering cells and genes...")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
logging.info(f"Shape after filtering: {adata.shape}")

# Verify counts layer survived filtering
logging.info("Verifying counts layer after filtering...")
verify_counts_layer(adata, "after filtering")

# Check if data is already normalized
logging.info("Checking normalization...")
totals = adata.X.sum(axis=1)
totals = totals.A1 if hasattr(totals, "A1") else np.array(totals).ravel()

if np.median(totals) > 2e4:
    logging.info("Normalizing data...")
    sc.pp.normalize_total(adata, target_sum=1e4)
else:
    logging.info("Data already normalized.")

# Check if data is log-transformed
logging.info("Checking log transformation...")
X_sample = adata.X[:100, :100]
arr = X_sample.A if hasattr(X_sample, "A") else X_sample

if np.allclose(arr, np.round(arr)):
    logging.info("Applying log transformation...")
    sc.pp.log1p(adata)
else:
    logging.info("Data already log-transformed.")

# Dimensionality reduction
logging.info("Calculating PCA and UMAP...")
sc.tl.pca(adata, n_comps=30, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
sc.tl.umap(adata)

# Save processed data
logging.info("Verifying counts layer before saving...")
verify_counts_layer(adata, "before saving")

adata.write_h5ad(processed_path)
logging.info(f"Saved processed adata to {processed_path}...")

# Verify the saved file has counts layer
logging.info("Verifying saved file has counts layer...")
adata_test = sc.read_h5ad(processed_path)
if "counts" in adata_test.layers:
    logging.info("✓ Counts layer successfully saved and verified!")
    logging.info(f"✓ Saved counts layer shape: {adata_test.layers['counts'].shape}")
else:
    logging.error("✗ Counts layer lost during save/load!")
    raise ValueError("Failed to save counts layer properly")
del adata_test  # Clean up memory

# Plot UMAP
logging.info("Plotting UMAPs...")
obs_lists = ["disease", "tissue_level_2", "cell_type"]
for obs in obs_lists:
    sc.pl.umap(
        adata,
        color=obs,
        title=f"{REF_DATASET} filtered and processed: UMAP colored by {obs}",
        save=f"umap_{REF_DATASET}_{obs}.pdf",
    )

logging.info("Filtering and processing script completed.")
