"""Image viewing module."""

# Import packages
import warnings
from logging import getLogger

import scanpy as sc
import seaborn as sns
import squidpy as sq

from recode_st.config import IOConfig, ViewImagesModuleConfig

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)


def run_view_images(config: ViewImagesModuleConfig, io_config: IOConfig):
    """Run the image viewing module."""
    # Set variables
    module_dir = io_config.output_dir / config.module_name

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set figure directory for this module (overrides global setting)
    sc.settings.figdir = module_dir

    # Set figure settings to ensure consistency across all modules
    # cmap = sns.color_palette("Spectral", as_cmap=True)
    cmap_blue = sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True)

    # Import data
    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(io_config.output_dir / "annotate" / "adata.h5ad")

    # View specific gene expression
    logger.info("Plotting genes of interest on tissue...")
    ROI_list = adata.obs["ROI"].unique().tolist()
    for roi in ROI_list:
        adata_roi = adata[adata.obs["ROI"] == roi]

        # Gene expression plots
        sq.pl.spatial_scatter(
            adata_roi,
            library_id="spatial",
            color=config.gene_list,
            cmap=cmap_blue,
            shape=None,
            size=1,
            img=False,
            save=f"gene_expression_{roi}.pdf",
        )
        logger.info(f"Saved gene expression plot for ROI {roi} to {module_dir}")

        # Cluster plots
        sq.pl.spatial_scatter(
            adata_roi,
            library_id="spatial",
            shape=None,
            outline=False,
            color=[
                "manual_annotation",
            ],  # TODO: remove hard coding
            wspace=0.4,
            size=1,
            save=f"clusters_{roi}.pdf",
            dpi=300,
        )
        logger.info(f"Saved cluster plot for ROI {roi} to {module_dir}")

    logger.info("Imaging module completed successfully.")
