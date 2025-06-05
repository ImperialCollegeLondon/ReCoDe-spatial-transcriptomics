"""Muspan module."""

# Import packages
import logging
import os
import warnings
from pathlib import Path

import scanpy as sc

from recode_st.paths import logging_path, output_path, zarr_path

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    # Set variables
    module_name = "6_muspan"  # name of the module
    module_dir = output_path / module_name

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set up logging
    os.makedirs(
        logging_path, exist_ok=True
    )  # should set up all these directories at the start of the pipeline?
    logging.basicConfig(
        filename=Path(logging_path) / f"{module_name}.txt",  # output file
        filemode="w",  # overwrites the file each time
        format="%(asctime)s - %(levelname)s - %(message)s",  # log format
        level=logging.INFO,  # minimum level to log
    )

    # change directory to output_path/module_name
    os.chdir(module_dir)
    logging.info(f"Changed directory to {module_dir}")

    # Import data
    logging.info("Loading Xenium data...")
    adata = sc.read_h5ad(output_path / "4_view_images/adata.h5ad")
    sdata = sd.read_zarr(zarr_path)
