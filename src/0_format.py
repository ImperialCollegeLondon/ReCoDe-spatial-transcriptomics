# Import packages
import os
import warnings  # ? what is the best way to suppress warnings from package inputs?
import logging
from pathlib import Path
import spatialdata as sd
from spatialdata_io import xenium

warnings.filterwarnings("ignore")

# Set directories
base_dir = Path(os.getenv("RECODE_BASE_DIR"))
input_path = base_dir
output_path = base_dir / "analysis"
logging_path = output_path / "logging"
xenium_path = input_path / "data/xenium"
zarr_path = base_dir / "data" / "xenium.zarr"

# Set up logging
logging.basicConfig(
    filename=logging_path / "0_format.txt",  # output file
    filemode="w",  # overwrites the file each time
    format="%(asctime)s - %(levelname)s - %(message)s",  # log format
    level=logging.INFO,  # minimum level to log
)

# Load into spatialdata format
logging.info("Reading Xenium data...")
sdata = xenium(xenium_path)

# Convert to zarr format
logging.info("Writing to Zarr...")
sdata.write(zarr_path)

# Convert to zarr format
logging.info("Finished formatting data.")
