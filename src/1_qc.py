"""Quality control module."""

# Import packages
import logging
import os
import warnings  # ? what is the best way to suppress warnings from package inputs?
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns
import spatialdata as sd

from .paths import base_dir, logging_path, output_path, zarr_path

warnings.filterwarnings("ignore")

# from helper_function import seed_everything # ? not working please help

# Set seed
# seed_everything(21122023) # ? not working please help

# Set variables
# ? How should I config this so a user can easily change them?
module_name = "1_qc"  # name of the module
min_counts = 10
min_cells = 5

# Confirm directories exist
if not Path(base_dir).exists():
    raise FileNotFoundError(f"Input path {base_dir} does not exist.")
if not Path(output_path).exists():
    raise FileNotFoundError(f"Output path {output_path} does not exist.")
if not Path(zarr_path).exists():
    raise FileNotFoundError(f"Zarr path {zarr_path} does not exist.")

# Create output directories if they do not exist
os.makedirs(Path(output_path) / module_name, exist_ok=True)

# Set up logging
os.makedirs(
    logging_path, exist_ok=True
)  # should set up all these directories at the start of the pipeline?
logging.basicConfig(
    filename=Path(logging_path) / f"{module_name}.txt",  # output file
    filemode="w",  # overwrites the file each time
    format="%(asctime)s - %(levelname)s - %(message)s",  # log format
    level=logging.INFO,  # minimum level to log
)

# Read in .zarr
logging.info("Loading Xenium data...")
sdata = sd.read_zarr(zarr_path)  #  read directly from the zarr store

logging.info("Done")

# # Save anndata object (stored in spatialdata.tables layer)
adata = sdata.tables["table"]  # contains the count matrix, cell and gene annotations

# $ Calculate and plot metrics

# Calculate quality control metrics
sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)

# Calculate percent negative DNA probe and percent negative decoding count
cprobes = (
    adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100
)
cwords = (
    adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100
)
logging.info(f"Negative DNA probe count % : {cprobes}")
logging.info(f"Negative decoding count % : {cwords}")

# Calculate averages
avg_total_counts = np.mean(adata.obs["total_counts"])
logging.info(f"Average number of transcripts per cell: {avg_total_counts}")

avg_total_unique_counts = np.mean(adata.obs["n_genes_by_counts"])
logging.info(f"Average unique transcripts per cell: {avg_total_unique_counts}")

area_max = np.max(adata.obs["cell_area"])
area_min = np.min(adata.obs["cell_area"])
logging.info(f"Max cell area: {area_max}")
logging.info(f"Min cell area: {area_min}")


# Plot
fig, axs = plt.subplots(1, 4, figsize=(15, 4))

axs[0].set_title("Total transcripts per cell")
sns.histplot(
    adata.obs["total_counts"],
    kde=False,
    ax=axs[0],
)

axs[1].set_title("Unique transcripts per cell")
sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    ax=axs[1],
)


axs[2].set_title("Area of segmented cells")
sns.histplot(
    adata.obs["cell_area"],
    kde=False,
    ax=axs[2],
)

axs[3].set_title("Nucleus ratio")
sns.histplot(
    adata.obs["nucleus_area"] / adata.obs["cell_area"],
    kde=False,
    ax=axs[3],
)

# Adjust layout and save the figure
plt.tight_layout()
plt.savefig(Path(output_path) / f"{module_name}/cell_summary_histograms.png", dpi=300)
plt.close()
logging.info(f"Saved plots to {module_name}/cell_summary_histograms.png")

# $ QC data #

# Filter cells
logging.info("Filtering cells and genes...")
sc.pp.filter_cells(adata, min_counts=min_counts)
sc.pp.filter_genes(adata, min_cells=min_cells)


# Normalize data
logging.info("Normalize data...")
adata.layers["counts"] = adata.X.copy()  # make copy of raw data
sc.pp.normalize_total(adata, inplace=True)  # normalize data
sc.pp.log1p(adata)  # Log transform data


# Save data
adata.write_h5ad(Path(output_path) / f"{module_name}/adata.h5ad")
logging.info(f"Data saved to {module_name}/adata.h5ad")
