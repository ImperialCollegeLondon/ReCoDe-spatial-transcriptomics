# Import packages
import warnings  # ? what is the best way to suppress warnings from package inputs?

warnings.filterwarnings("ignore")

import numpy as np

import spatialdata as sd
from spatialdata_io import xenium

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq

from pathlib import Path

import torch
import os

import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler("my_script.log"), logging.StreamHandler()],
)


# from helper_function import seed_everything # ? not working please help

# Set seed
# seed_everything(21122023) # ? not working please help

# Set variables
# ? How should I config this so a user can easily change them?
min_counts = 10
min_cells = 5


# Set directories -
input_path = "/Users/sarapatti/Desktop/PhD_projects/Llyod_lab/ReCoDe-spatial-transcriptomics"
output_path = "/Users/sarapatti/Desktop/PhD_projects/Llyod_lab/ReCoDe-spatial-transcriptomics/analysis"
xenium_path = Path(input_path) / "data/xenium"
zarr_path = Path(input_path) / "data/xenium.zarr"

# Load into spatialdata format
logging.info("Reading Xenium data...")
sdata = xenium(xenium_path)

# Convert to zarr format
logging.info("Writing to Zarr...")
sdata.write(zarr_path)

# Read in .zarr
logging.info("Loading Xenium data...")
sdata = sd.read_zarr(zarr_path)  #  read directly from the zarr store

logging.info("Done")

# # Save anndata oject (stored in spatialdata.tables layer)
adata = sdata.tables[
    "table"
]  # contains the count matrix, cell and gene annotations

# $ Calculate and plot metrics

# Calculate quality control metrics
sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)

# Calculate percent negative DNA probe and percent negative decoding count
cprobes = (
    adata.obs["control_probe_counts"].sum()
    / adata.obs["total_counts"].sum()
    * 100
)
cwords = (
    adata.obs["control_codeword_counts"].sum()
    / adata.obs["total_counts"].sum()
    * 100
)
print(f"Negative DNA probe count % : {cprobes}")
print(f"Negative decoding count % : {cwords}")

# Calculate averages
avg_total_counts = np.mean(adata.obs["total_counts"])
print(f"Average number of transcripts per cell: {avg_total_counts}")

avg_total_unique_counts = np.mean(adata.obs["n_genes_by_counts"])
print(f"Average unique transcripts per cell: {avg_total_unique_counts}")

area_max = np.max(adata.obs["cell_area"])
area_min = np.min(adata.obs["cell_area"])
print(f"Max cell area: {area_max}")
print(f"Min cell area: {area_min}")


# Plot
# ^ How do I save these plots?
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

# $ QC data #

# Filter cells
sc.pp.filter_cells(adata, min_counts=min_counts)
sc.pp.filter_genes(adata, min_cells=min_cells)


# Normalize data
adata.layers["counts"] = adata.X.copy()  # make copy of raw data
sc.pp.normalize_total(adata, inplace=True)  # normalize data
# Log transform data
sc.pp.log1p(adata)


# Save data
# ? Should I have the QC data as a few .5had file each time?
