"""Image viewing module."""

# Import packages
import warnings
from logging import getLogger

import scanpy as sc
import seaborn as sns
import squidpy as sq

from recode_st.config import IOConfig, ViewImagesModuleConfig
from recode_st.helper_function import configure_scanpy_figures

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)


def run_view_images(config: ViewImagesModuleConfig, io_config: IOConfig):
    """Run the image viewing module."""
    # Set variables
    module_dir = io_config.output_dir / config.module_name

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set figure settings to ensure consistency across all modules
    configure_scanpy_figures(str(io_config.output_dir))
    cmap = sns.color_palette("Spectral", as_cmap=True)

    # Import data
    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(io_config.output_dir / "3_annotate" / "adata.h5ad")

    # View plots
    logger.info("Visualize clusters on tissue...")
    sq.pl.spatial_scatter(
        adata,
        library_id="spatial",
        shape=None,
        outline=False,
        color=["leiden", "total_counts"],
        wspace=0.4,
        size=1,
        save=module_dir / "leiden_clusters.png",
        dpi=300,
    )
    logger.info(f"Saved leiden clusters plot to {module_dir / 'leiden_clusters.png'}")

    # View specific gene expression
    logger.info("Plotting genes of interest on tissue...")
    sq.pl.spatial_scatter(
        adata,
        library_id="spatial",
        color=config.gene_list,
        shape=None,
        size=2,
        img=False,
        save=module_dir / "gene_expression.png",
    )
    logger.info(f"Saved gene expression plot to {module_dir / 'gene_expression.png'}")

    # Save anndata object
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")
    logger.info("Imaging module completed successfully.")
