"""Image viewing module."""

# Import packages
import os
import warnings
from logging import getLogger

import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

from recode_st.logging_config import configure_logging
from recode_st.paths import output_path

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def run_view_images():
    """Run the image viewing module."""
    # Set variables
    module_name = "4_view_images"  # name of the module
    gene_list = ["EPCAM", "CD3D", "CD68", "VWF", "PTPRC", "ACTA2"]
    module_dir = output_path / module_name

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # change directory to output_path/module_name
    os.chdir(module_dir)
    logger.info(f"Changed directory to {module_dir}")

    # Import data
    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(module_dir / "3_annotate/adata.h5ad")

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
    )
    plt.tight_layout()
    plt.savefig(module_dir / f"{module_name}/images.png", dpi=300)
    plt.close()
    logger.info(f"Saved plots to {module_dir / 'images.png'}")

    # View specific gene expression
    logger.info("Plotting genes of interest on tissue...")
    sq.pl.spatial_scatter(
        adata,
        library_id="spatial",
        color=gene_list,
        shape=None,
        size=2,
        img=False,
        save="gene_expression.png",
    )

    # Save anndata object
    adata.write_h5ad(module_dir / f"{module_name}/adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")
    logger.info("Imaging module completed successfully.")


if __name__ == "__main__":
    # Set up logger
    configure_logging()
    logger = getLogger("recode_st.4_view_images")

    try:
        run_view_images()
    except FileNotFoundError as err:
        logger.error(f"File not found: {err}")
