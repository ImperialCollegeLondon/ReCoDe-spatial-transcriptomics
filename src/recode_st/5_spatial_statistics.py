"""Spatial statistics module."""

import os
import warnings
from logging import getLogger
from pathlib import Path

import scanpy as sc
import squidpy as sq

from recode_st.logging_config import configure_logging
from recode_st.paths import output_path

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def run_spatial_statistics():
    """Run spatial statistics on Xenium data."""
    # Set variables
    module_name = "5_spatial_stats"  # name of the module
    module_dir = output_path / module_name

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # change directory to output_path/module_name
    os.chdir(module_dir)
    logger.info(f"Changed directory to {module_dir}")

    # Import data
    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(module_dir / "4_view_images/adata.h5ad")

    # $ Calculate spatial statistics

    logger.info("Building spatial neighborhood graph...")
    sq.gr.spatial_neighbors(
        adata, coord_type="generic", delaunay=True
    )  # compute connectivity

    logger.info("Computing and plotting centrality scores...")
    sq.gr.centrality_scores(adata, cluster_key="leiden")
    sq.pl.centrality_scores(
        adata, cluster_key="leiden", figsize=(16, 5), save="_plot.png"
    )

    # # $ Compute co-occurrence probability
    # logger.info("Computing co-occurrence probability...")
    # # Create subset table layer
    # adata_subsample = sc.pp.subsample(
    #     adata, fraction=0.5, copy=True
    # )  # subsample to speed up computation

    # # Visualize co-occurrence
    # sq.gr.co_occurrence(
    #     adata_subsample,
    #     cluster_key="leiden",
    # )
    # sq.pl.co_occurrence(
    #     adata_subsample,
    #     cluster_key="leiden",
    #     clusters="12",
    #     figsize=(10, 10),
    #     save="_plot.png",
    # )

    # $ Neighborhood enrichment analysis
    logger.info("Performing neighborhood enrichment analysis...")
    sq.gr.nhood_enrichment(adata, cluster_key="leiden")

    # Plot neighborhood enrichment
    sq.pl.nhood_enrichment(
        adata,
        cluster_key="leiden",
        figsize=(8, 8),
        title="Neighborhood enrichment adata",
        save="_plot.png",
    )

    # $ Moran's I
    logger.info("Calculating Moran's I...")

    # # Build spatial neighborhood graph on a subsampled dataset
    # sq.gr.spatial_neighbors(adata_subsample, coord_type="generic", delaunay=True)

    # # Calculate Moran's I for spatial autocorrelation on subsampled data
    # sq.gr.spatial_autocorr(
    #     adata_subsample,
    #     mode="moran",
    #     n_perms=100,
    #     n_jobs=1,
    # )

    # # Save Moran's I results
    # adata_subsample.uns["moranI"].to_csv(
    #     Path(output_path) / module_name / "moranI_results.csv",
    #     index=True,
    # )

    # Save anndata object
    adata.write_h5ad(Path(output_path) / f"{module_name}/adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")
    logger.info("Spatial statistics module completed successfully.")


if __name__ == "__main__":
    # Set up logger
    configure_logging()
    logger = getLogger("recode_st.5_spatial_statistics")

    try:
        run_spatial_statistics()
    except FileNotFoundError as err:
        logger.error(f"File not found: {err}")
