"""Module for defining inputs and outputs paths."""

from pathlib import Path

# Set directories
base_dir = Path(__file__).parent.parent.parent  # TODO: Make this configurable
data_path = base_dir / "data"
xenium_path = data_path / "xenium"
zarr_path = data_path / "xenium.zarr"
output_path = base_dir / "analysis"
logging_path = output_path / "logs"

# Confirm directories exist
if not Path(base_dir).exists():
    raise FileNotFoundError(f"Input path {base_dir} does not exist.")
if not Path(data_path).exists():
    raise FileNotFoundError(f"Input path {base_dir} does not exist.")
if not Path(zarr_path).exists():
    raise FileNotFoundError(f"Input path {base_dir} does not exist.")

# Make sure output directories exist
if not Path(output_path).exists():
    output_path.mkdir(parents=True, exist_ok=True)
if not Path(logging_path).exists():
    logging_path.mkdir(parents=True, exist_ok=True)
