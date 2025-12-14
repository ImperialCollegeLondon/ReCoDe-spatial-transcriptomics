"""Spatial statistics module."""

import warnings
from logging import getLogger

import scanpy as sc
import seaborn as sns
import squidpy as sq

from recode_st.config import IOConfig, SpatialStatisticsModuleConfig
from recode_st.helper_function import configure_scanpy_figures

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)


def run_spatial_statistics(config: SpatialStatisticsModuleConfig, io_config: IOConfig):
    """Run spatial statistics on Xenium data."""
    # Set variables
    module_dir = io_config.output_dir / config.module_name

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set figure settings to ensure consistency across all modules
    configure_scanpy_figures(str(io_config.output_dir))
    cmap = sns.color_palette("Spectral", as_cmap=True)

    # Import data
    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(io_config.output_dir / "4_view_images" / "adata.h5ad")

    # Calculate spatial statistics
    logger.info("Building spatial neighborhood graph...")
    sq.gr.spatial_neighbors(
        adata, coord_type="generic", delaunay=True
    )  # compute connectivity

    logger.info("Computing and plotting centrality scores...")
    sq.gr.centrality_scores(adata, cluster_key="leiden")
    sq.pl.centrality_scores(
        adata,
        cluster_key="leiden",
        figsize=(16, 5),
        save=module_dir / "centrality_scores.png",
    )
    logger.info(
        f"Centrality scores plot saved to {module_dir / 'centrality_scores.png'}"
    )

    # Compute co-occurrence probability
    logger.info("Computing co-occurrence probability...")
    # Create subset table layer
    adata_subsample = sc.pp.subsample(
        adata, fraction=0.5, copy=True
    )  # subsample to speed up computation

    # Visualize co-occurrence
    sq.gr.co_occurrence(
        adata_subsample,
        cluster_key="leiden",
    )
    sq.pl.co_occurrence(
        adata_subsample,
        cluster_key="leiden",
        clusters="12",
        figsize=(10, 10),
        save=module_dir / "co_occurrence.png",
    )
    logger.info(f"Co-occurrence plot saved to {module_dir / 'co_occurrence.png'}")

    # Neighborhood enrichment analysis
    logger.info("Performing neighborhood enrichment analysis...")
    sq.gr.nhood_enrichment(adata, cluster_key="leiden")

    # Plot neighborhood enrichment
    sq.pl.nhood_enrichment(
        adata,
        cluster_key="leiden",
        figsize=(8, 8),
        title="Neighborhood enrichment adata",
        save=module_dir / "nhood_enrichment.png",
    )
    logger.info(
        f"Neighborhood enrichment plot saved to {module_dir / 'nhood_enrichment.png'}"
    )

    # Moran's I
    logger.info("Calculating Moran's I...")

    # Build spatial neighborhood graph on a subsample dataset
    sq.gr.spatial_neighbors(adata_subsample, coord_type="generic", delaunay=True)

    # Calculate Moran's I for spatial autocorrelation on subsample data
    sq.gr.spatial_autocorr(
        adata_subsample,
        mode="moran",
        n_perms=100,
        n_jobs=1,
    )

    # Save Moran's I results
    adata_subsample.uns["moranI"].to_csv(module_dir / "moranI_results.csv", index=True)
    logger.info(f"Moran's I results saved to {module_dir / 'moranI_results.csv'}")

    # Save anndata object
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")
    logger.info("Spatial statistics module completed successfully.")
