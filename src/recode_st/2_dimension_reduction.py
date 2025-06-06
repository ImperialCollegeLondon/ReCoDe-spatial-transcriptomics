"""Dimension reduction module."""

import os
import warnings
from logging import getLogger

import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

from recode_st.helper_function import seed_everything
from recode_st.logging_config import configure_logging
from recode_st.paths import output_path

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def run_dimension_reduction():
    """Run dimension reduction on Xenium data."""
    # Set variables
    module_name = "2_DR"
    module_dir = output_path / module_name

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Import data
    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(output_path / "1_qc/adata.h5ad")

    # Perform dimension reduction analysis
    logger.info("Compute PCA...")
    sc.pp.pca(adata)  # compute principal components
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)
    plt.tight_layout()
    plt.savefig(
        module_dir / "pca_variance_ratio.png"
    )  # save the figure with the module name

    logger.info("Compute neighbors...")
    sc.pp.neighbors(adata)  # compute a neighborhood graph

    logger.info("Xreate UMAPs and cluster cells..")
    sc.tl.umap(adata)  # calculate umap
    sc.tl.leiden(
        adata,
        resolution=1.0,  # choose resolution for clustering
        key_added="leiden",
    )  # name leiden clusters

    # change directory to output_path/module_name
    os.chdir(module_dir)
    logger.info(f"Changed directory to {module_dir}")

    # plot UMAP
    logger.info("Plotting UMAPs...")
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
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")


if __name__ == "__main__":
    # Set up logger
    configure_logging()
    logger = getLogger("recode_st.2_dimension_reduction")

    # Set seed
    seed_everything(21122023)

    try:
        run_dimension_reduction()
    except FileNotFoundError as err:
        logger.error(f"File not found: {err}")
