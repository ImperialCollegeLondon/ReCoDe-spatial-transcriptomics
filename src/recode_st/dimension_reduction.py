"""Dimension reduction module."""

import warnings
from logging import getLogger
from pathlib import Path
from typing import Any, Literal

import geosketch
import numpy as np
import scanpy as sc
import seaborn as sns
import squidpy as sq

from recode_st.config import DimensionReductionModuleConfig, IOConfig
from recode_st.helper_function import configure_scanpy_figures

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)

SubsampleStrategy = Literal["none", "compute", "load"]


def create_subsample(
    adata: sc.AnnData,
    n_total: int = 10_000,
    min_cells_per_roi: int = 50,
    n_pca: int = 50,
) -> sc.AnnData:
    """Create a subsampled dataset using geosketch.

    Args:
        adata: Input AnnData object (full dataset)
        n_total: Total number of cells to subsample
        min_cells_per_roi: Minimum cells to keep per ROI
        n_pca: Number of PCA components for sketching

    Returns:
        Subsampled AnnData object
    """
    logger.info("Subsampling data...")
    orig_size = len(adata)
    logger.info(f"Original size: {orig_size} cells")

    # Ensure PCA is computed
    if "X_pca" not in adata.obsm:
        logger.info("Computing PCA for sketching...")
        sc.pp.pca(adata, n_comps=n_pca)

    roi_counts = adata.obs["ROI"].value_counts()
    sampled_indices = []
    logger.info(f"Sketching {n_total} cells across {len(roi_counts)} ROIs...")

    for roi in adata.obs["ROI"].unique():
        roi_mask = adata.obs["ROI"] == roi
        roi_data = adata[roi_mask]

        # Calculate proportional sample size for this ROI
        n_roi_sketch = int(n_total * (roi_counts[roi] / len(adata)))
        n_roi_sketch = max(min_cells_per_roi, n_roi_sketch)
        n_roi_sketch = min(n_roi_sketch, len(roi_data))

        logger.info(f"  ROI {roi}: sketching {n_roi_sketch} from {len(roi_data)} cells")

        # Sketch or take all
        if n_roi_sketch >= len(roi_data):
            roi_sketch_local = np.arange(len(roi_data))
        else:
            roi_sketch_local = geosketch.gs(
                roi_data.obsm["X_pca"], n_roi_sketch, replace=False
            )

        global_indices = np.where(roi_mask)[0][roi_sketch_local]
        sampled_indices.extend(global_indices)

    adata_sub = adata[sampled_indices, :].copy()

    # Log results
    logger.info("\nSubsampling Results:")
    logger.info(f"  Original: {orig_size} cells")
    logger.info(f"  Subsampled: {len(adata_sub)} cells")
    logger.info(f"  Reduction: {100 * (1 - len(adata_sub) / orig_size):.1f}%")
    logger.info("\nCells per ROI:")
    logger.info(adata_sub.obs["ROI"].value_counts().sort_index())

    return adata_sub


def load_data(
    io_config: IOConfig,
    norm_approach: str,
    subsample_strategy: SubsampleStrategy,
    subsample_path: Path,
    n_total: int = 10_000,
    min_cells_per_roi: int = 50,
    n_pca: int = 50,
) -> sc.AnnData:
    """Load data based on subsampling strategy.

    Args:
        io_config: IO configuration
        norm_approach: Normalization approach label
        subsample_strategy: One of "none", "compute", or "load"
        subsample_path: Path to save/load subsampled data
        n_total: Total cells for subsampling
        min_cells_per_roi: Minimum cells per ROI
        n_pca: Number of PCA components

    Returns:
        AnnData object ready for analysis

    Raises:
        FileNotFoundError: If load strategy is used but file doesn't exist
        ValueError: If invalid strategy is provided
    """
    data_path = (
        io_config.output_dir / "1_quality_control" / f"adata_{norm_approach}.h5ad"
    )

    if subsample_strategy == "none":
        logger.info("Loading full dataset (no subsampling)...")
        adata = sc.read_h5ad(data_path)

        # Compute PCA if needed
        if "X_pca" not in adata.obsm:
            logger.info("Computing PCA...")
            sc.pp.pca(adata, n_comps=n_pca)

        logger.info(f"Loaded {len(adata)} cells")
        return adata

    elif subsample_strategy == "compute":
        logger.info("Creating new subsample...")

        # Load full dataset
        adata = sc.read_h5ad(data_path)

        # Create subsample
        adata_sub = create_subsample(
            adata, n_total=n_total, min_cells_per_roi=min_cells_per_roi, n_pca=n_pca
        )

        # Save subsample
        logger.info(f"Saving subsample to {subsample_path}")
        adata_sub.write_h5ad(subsample_path)

        return adata_sub

    elif subsample_strategy == "load":
        logger.info("Loading existing subsample...")

        if not subsample_path.exists():
            raise FileNotFoundError(
                f"Subsample file not found: {subsample_path}\n"
                f"Run with subsample_strategy='compute' first."
            )

        adata = sc.read_h5ad(subsample_path)
        logger.info(f"Loaded {len(adata)} cells")
        logger.info("Cells per ROI:")
        logger.info(adata.obs["ROI"].value_counts().sort_index())

        return adata

    else:
        raise ValueError(
            f"Invalid subsample_strategy: '{subsample_strategy}'. "
            f"Must be one of: 'none', 'compute', 'load'"
        )


def compute_dimensionality_reduction(
    adata: sc.AnnData,
    n_pca: int = 30,  # number of PCA components to use to compute neighbors
    n_neighbors: int = 10,  # number of neighbors for graph
    min_dist: float = 0.5,  # minimum distance for UMAP
    resolution: float = 1,  # resolution for Leiden clustering
    cluster_name: str = "leiden",  # name for cluster annotation
) -> sc.AnnData:
    """Compute PCA, neighbors, UMAP, and clustering.

    Args:
        adata: Input AnnData object
        n_pca: Number of PCA components to use
        n_neighbors: Number of neighbors for graph
        resolution: Resolution for Leiden clustering
        cluster_name: Name for cluster annotation
        min_dist: Minimum distance for UMAP

    Returns:
        Nothing. Saves figures and AnnData with computed dimensionality reduction.
    """
    logger.info(f"Computing neighbors (k={n_neighbors})...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pca)

    logger.info("Computing UMAP...")
    sc.tl.umap(adata, min_dist=min_dist)

    logger.info(f"Clustering with resolution={resolution}...")
    sc.tl.leiden(adata, resolution=resolution, key_added=cluster_name)

    return adata


def plot_dimensionality_reduction(
    adata: sc.AnnData,
    norm_approach: str,
    module_name: str,
    n_neighbors: int,
    figdir: Path,
    cmap: Any = sns.color_palette("crest", as_cmap=True),
    cluster_name: str = "leiden",
    config: DimensionReductionModuleConfig = None,
) -> None:
    """Create and save dimensionality reduction plots.

    Args:
        adata: AnnData with computed UMAP and clusters
        norm_approach: Normalization approach label
        module_name: Name of module for file naming
        n_neighbors: Number of neighbors used
        figdir: Directory to save figures
        cmap: Colormap for plots
        cluster_name: Name of cluster annotation
        config: Configuration object containing visualization settings
    """
    logger.info("Plotting UMAPs...")

    # QC and metadata plots
    logger.info("Plotting UMAPs with QC meterics...")
    sc.pl.umap(
        adata,
        color=["total_counts", "n_genes_by_counts", cluster_name],
        ncols=3,
        cmap=cmap,
        wspace=0.4,
        show=False,
        save=f"_{module_name}_{norm_approach}_neighbors_{n_neighbors}.png",
        frameon=False,
    )

    logger.info("Plotting UMAPs with observation fields...")
    sc.pl.umap(
        adata,
        color=config.obs_vis_list,
        ncols=3,
        cmap=cmap,
        wspace=0.4,
        show=False,
        save=f"_{module_name}_{norm_approach}_neighbors_{n_neighbors}.png",
        frameon=False,
    )

    logger.info("Plotting UMAPs with marker genes...")
    sc.pl.umap(
        adata,
        color=config.marker_genes,
        cmap=cmap,
        ncols=4,
        wspace=0.4,
        show=False,
        save=f"_{module_name}_cell_markers_{norm_approach}_neighbors_{n_neighbors}.png",
        frameon=False,
    )

    logger.info(f"UMAP plots saved to {figdir}")


def plot_spatial_distribution(
    adata: sc.AnnData,
    cluster_name: str,
    module_dir: Path,
) -> None:
    """Plot spatial distribution of clusters for each ROI.

    Args:
        adata: AnnData with cluster annotations
        cluster_name: Name of cluster annotation
        module_dir: Directory to save plots
    """
    logger.info(f"Plotting {cluster_name} spatial distribution...")

    for roi in adata.obs["ROI"].unique():
        subset = adata[adata.obs["ROI"] == roi]
        sq.pl.spatial_scatter(
            subset,
            library_id="spatial",
            shape=None,
            color=[cluster_name],
            wspace=0.4,
            save=module_dir / f"{cluster_name}_{roi}_spatial.png",
        )

    logger.info(f"Spatial plots saved to {module_dir}")


def run_dimension_reduction(
    config: DimensionReductionModuleConfig, io_config: IOConfig
) -> sc.AnnData:
    """Run dimension reduction on Xenium data.

    Args:
        config: Configuration for dimension reduction
        io_config: IO configuration

    Returns:
        AnnData object with computed dimensionality reduction
    """
    # Setup
    module_dir = io_config.output_dir / config.module_name
    module_dir.mkdir(exist_ok=True, parents=True)

    # Set figure directory for this module (overrides global setting)
    sc.settings.figdir = module_dir

    # Define variables
    CLUSTER_NAME = "leiden"

    # Set figure settings to ensure consistency across all modules
    configure_scanpy_figures(str(io_config.output_dir))
    cmap = sns.color_palette("Spectral", as_cmap=True)

    # Load data based on subsampling strategy
    subsample_path = module_dir / f"adata_subsample_{config.norm_approach}.h5ad"

    adata = load_data(
        io_config=io_config,
        norm_approach=config.norm_approach,
        subsample_strategy=config.subsample_strategy,
        subsample_path=subsample_path,
        n_total=config.subsample_n_total,
        min_cells_per_roi=config.subsample_min_cells_per_roi,
        n_pca=config.n_pca,
    )

    # Plot PCA variance
    logger.info("Plotting PCA variance...")
    sc.pl.pca_variance_ratio(
        adata,
        log=True,
        n_pcs=80,
        show=False,
        save=f"_{config.module_name}.png",
    )

    # Compute dimensionality reduction
    logger.info("Computing dimension reduction...")
    adata = compute_dimensionality_reduction(
        adata=adata,
        n_pca=config.n_pca,
        n_neighbors=config.n_neighbors,
        resolution=config.resolution,
        cluster_name=CLUSTER_NAME,
    )

    # Save final results
    output_path = module_dir / "adata.h5ad"
    adata.write_h5ad(output_path)
    logger.info(f"\nFinal results saved to {output_path}")

    logger.info("Plotting dimension reduction...")

    # Plot results
    plot_dimensionality_reduction(
        adata=adata,
        cluster_name=CLUSTER_NAME,
        norm_approach=config.norm_approach,
        module_name=config.module_name,
        n_neighbors=config.n_neighbors,
        cmap=cmap,
        figdir=module_dir,
    )

    logger.info("Plotting spatial distribution of clusters...")
    # Plot spatial distribution of clusters
    spatial_plots_dir = module_dir / "spatial_plots"
    spatial_plots_dir.mkdir(exist_ok=True, parents=True)
    plot_spatial_distribution(
        adata=adata,
        cluster_name=CLUSTER_NAME,
        module_dir=spatial_plots_dir,
    )

    logger.info(f"\nDimension reduction module '{config.module_name}' complete.\n")
