"""Muspan module."""

import os
import warnings
from logging import getLogger

import scanpy as sc
import spatialdata as sd

from recode_st.logging_config import configure_logging
from recode_st.paths import output_path, zarr_path

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def run_muspan():
    """Run Muspan on Xenium data."""
    # Set variables
    module_name = "6_muspan"  # name of the module
    module_dir = output_path / module_name

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # change directory to output_path/module_name
    os.chdir(module_dir)
    logger.info(f"Changed directory to {module_dir}")

    # Import data
    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(output_path / "4_view_images/adata.h5ad")  # noqa: F841
    sdata = sd.read_zarr(zarr_path)  # noqa: F841


if __name__ == "__main__":
    # Set up logger
    configure_logging()
    logger = getLogger("recode_st.6_muspan")  # re-name the logger to match the module

    try:
        run_muspan()
    except FileNotFoundError as err:
        logger.error(f"File not found: {err}")
