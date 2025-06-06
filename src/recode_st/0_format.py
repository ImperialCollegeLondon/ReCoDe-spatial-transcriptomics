"""Module for formatting Xenium data into Zarr format."""

import warnings
from logging import getLogger
from pathlib import Path

from spatialdata_io import xenium

from recode_st.logging_config import configure_logging
from recode_st.paths import xenium_path, zarr_path

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def convert_xenium_to_zarr(xenium_path: Path, zarr_path: Path):
    """Convert Xenium data to Zarr format."""
    # Load Xenium data
    logger.info("Reading Xenium data...")
    sdata = xenium(xenium_path)

    # Write to Zarr format
    logger.info("Writing to Zarr...")
    sdata.write(zarr_path)

    logger.info("Finished formatting data.")


if __name__ == "__main__":
    # Set up logger
    configure_logging()
    logger = getLogger("recode_st.0_format")  # re-name the logger to match the module

    try:
        convert_xenium_to_zarr(xenium_path, zarr_path)
    except FileNotFoundError as err:
        logger.error(f"File not found: {err}")
