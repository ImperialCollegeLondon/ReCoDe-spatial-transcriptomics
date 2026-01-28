"""Import and format Xenium data into Zarr and AnnData formats."""

import warnings
from logging import getLogger
from pathlib import Path

import anndata as ad
from spatialdata_io import xenium

from recode_st.config import IOConfig

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)


def extract_roi_name(folder_name: str) -> str:
    """Extract ROI name from folder name."""
    # Example: "output-...__IPF_RBH_19__20251001__141533" → "IPF_RBH_19"
    parts = folder_name.split("__")
    for part in parts:
        if part.startswith(("IPF", "COPD", "PM08")):
            return part
    raise ValueError(f"Could not parse ROI name from: {folder_name}")


def convert_all_xenium(io_config: IOConfig):
    """Convert all Xenium ROIs into Zarrs and one combined AnnData."""
    xenium_path = Path(io_config.xenium_dir)

    zarr_root = Path(io_config.zarr_dir)
    zarr_root.mkdir(parents=True, exist_ok=True)

    output_data_dir = Path(io_config.output_data_dir)
    output_data_dir.mkdir(parents=True, exist_ok=True)

    adata_dir = Path(io_config.adata_dir)
    adata_dir.mkdir(parents=True, exist_ok=True)

    all_adatas = []

    print(f"Starting conversion of Xenium data in {xenium_path}")

    if not xenium_path.exists():
        raise FileNotFoundError(f"Xenium root directory not found: {xenium_path}")

    # Loop over date/run folders directly inside Xenium root
    for date_dir in xenium_path.iterdir():
        if not date_dir.is_dir():
            continue

        run_name = date_dir.name  # e.g. "20251001__141239__SP25164_SARA_PATTI_RUN_1"
        logger.info(f"Processing run: {run_name}")

        # Loop over ROI folders
        for roi_folder in date_dir.iterdir():
            if roi_folder.is_dir() and roi_folder.name.startswith("output-"):
                roi_name = extract_roi_name(roi_folder.name)
                logger.info(f"Processing ROI: {roi_name}")

                try:
                    sdata = xenium(roi_folder)
                    logger.info(f"Loaded Xenium data for {roi_name}")
                except Exception as err:
                    logger.error(f"Failed loading Xenium data for {roi_name}: {err}")
                    continue

                # Save Zarr per ROI
                roi_zarr_path = zarr_root / f"{roi_name}.zarr"
                try:
                    sdata.write(roi_zarr_path, overwrite=True)
                    logger.info(f"Zarr for {roi_name} saved to {roi_zarr_path}")
                except Exception as err:
                    logger.error(f"Failed writing Zarr for {roi_name}: {err}")
                    continue

                # Extract AnnData and add labels
                try:
                    # extract AnnData from SpatialData
                    adata = sdata.tables["table"]

                    # Add metadata to AnnData
                    adata.obs["ROI"] = roi_name  # add ROI name
                    adata.obs["batch"] = run_name  # add run/date name
                    condition = (  # Add condition "IPF", "COPD", "PM08", "Unknown"
                        roi_name.split("_")[0]
                        if roi_name.split("_")[0] in ["IPF", "COPD", "PM08"]
                        else "Unknown"
                    )
                    adata.obs["condition"] = condition
                    logger.info(
                        f"Added {roi_name}, {run_name}, and {condition} "
                        f"to AnnData for {roi_name}"
                    )
                    # Rename obs_names to include ROI so they are unique across ROIs
                    adata.obs_names = roi_name + "_" + adata.obs_names.astype(str)

                    if adata.obs_names.is_unique is False:
                        raise ValueError(
                            f"Non-unique obs_names in AnnData for {roi_name}"
                        )
                    else:
                        logger.info(f"obs_names are unique in AnnData for {roi_name}")
                        logger.info(f"obs_names: {adata.obs_names[:3]} ...")
                        logger.info(
                            f"AnnData for {roi_name} has shape {adata.shape} "
                            f"with {adata.n_vars} genes."
                        )

                    # Save individual AnnData
                    adata.write(adata_dir / f"{roi_name}.h5ad")
                    logger.info(f"AnnData for {roi_name} saved.")
                    all_adatas.append(adata)
                except KeyError:
                    logger.warning(f"No table found in {roi_name}, skipping AnnData.")

    # Concatenate all adatas
    if all_adatas:
        combined = ad.concat(all_adatas, join="outer", fill_value=0)
        combined_path = adata_dir / "all_samples.h5ad"
        combined.write(combined_path)
        logger.info(f"Combined AnnData written to {combined_path}")
    else:
        logger.warning("No AnnData objects were created. Combined AnnData not written.")


def run_format(io_config: IOConfig):
    """Run the formatting module across multiple ROIs."""
    convert_all_xenium(io_config)
    logger.info("Finished formatting all ROIs.")
