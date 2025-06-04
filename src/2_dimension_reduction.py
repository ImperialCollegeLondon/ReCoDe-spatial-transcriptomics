# Import packages
import warnings  # ? what is the best way to suppress warnings from package inputs?

import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq
import torch

warnings.filterwarnings("ignore")

# from helper_function.py import seed_everything

# Set seed
# seed_everything(21122023)

# Set variables
module_name = "2_DR"

# Set directories
base_dir = Path(os.getenv("RECODE_BASE_DIR"))
input_path = base_dir
output_path = base_dir / "analysis"
logging_path = output_path / "logging"

# Confirm directories exist
if not Path(input_path).exists():
    raise FileNotFoundError(f"Input path {input_path} does not exist.")
if not Path(output_path).exists():
    raise FileNotFoundError(f"Output path {output_path} does not exist.")

# Create output directories if they do not exist
os.makedirs(output_path / module_name, exist_ok=True)

# Set up logging
os.makedirs(logging_path, exist_ok=True)
logging.basicConfig(
    filename=logging_path / f"{module_name}.txt",  # output file
    filemode="w",  # overwrites the file each time
    format="%(asctime)s - %(levelname)s - %(message)s",  # log format
    level=logging.INFO,  # minimum level to log
)

# Import data
logging.info("Loading Xenium data...")
adata = sc.read_h5ad(output_path / "1_qc/adata.h5ad")

# Preform dimension reduction analysis
logging.info("Compute PCA...")
sc.pp.pca(adata)  # compute principal components
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)
plt.tight_layout()
plt.savefig(
    output_path / module_name / "pca_variance_ratio.png"
)  # save the figure with the module name

logging.info("Compute neighbors...")
sc.pp.neighbors(adata)  # compute a neighborhood graph

logging.info("Xreate UMAPs and cluster cells..")
sc.tl.umap(adata)  # calculate umap
sc.tl.leiden(
    adata,
    resolution=1.0,  # choose resolution for clustering
    key_added="leiden",
)  # name leiden clusters


# change directory to output_path/module_name
os.chdir(output_path / module_name)

# plot UMAP
logging.info("Plotting UMAPs...")
sc.pl.umap(
    adata,
    color=[
        "total_counts",
        "n_genes_by_counts",
        "leiden",
    ],
    wspace=0.4,
    show=False,
    save=f"_{module_name}.png",  # save the figure with the module name
    frameon=False,
)

# plot visualization of leiden clusters
sq.pl.spatial_scatter(
    adata,
    library_id="spatial",
    shape=None,
    color=[
        "leiden.0",
    ],
    wspace=0.4,
    save="leiden_spatial.png",
)

# Save anndata object
adata.write_h5ad(output_path / f"{module_name}/adata.h5ad")
logging.info(f"Data saved to {module_name}/adata.h5ad")
