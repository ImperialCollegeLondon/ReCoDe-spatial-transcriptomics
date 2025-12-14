"""Dimension reduction module."""

import warnings
from logging import getLogger

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from zarr.errors import PathNotFoundError

from recode_st.config import Config, IOConfig, SubsampleModuleConfig
from recode_st.helper_function import configure_scanpy_figures

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)


def run_subsampled_data(
    config: SubsampleModuleConfig, io_config: IOConfig, main_config: Config
):
    """Run data subsampling on Xenium data."""
    # Set variables
    module_dir = io_config.output_dir / config.module_name

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set figure directory for this module (overrides global setting)
    sc.settings.figdir = module_dir

    # Ensure global visualization settings are applied
    # This ensures consistency across all modules
    configure_scanpy_figures(str(io_config.output_dir))

    # Import data
    try:
        logger.info("Loading Xenium data...")
        combined_path = io_config.adata_dir / "all_samples.h5ad"

        # Read the file normally - we need the data in memory for QC calculations
        adata = sc.read_h5ad(combined_path)

    except PathNotFoundError as err:
        logger.error(f"File not found (or not a valid AnnData file): {combined_path}")
        raise err

    logger.info("Xenium data loaded.")

    logger.info(f"Dataset contains {len(adata)} cells prior to subsetting.")

    logger.info("Subsampling data by ROI...")
    # Subsample data
    adata = subsample_by_roi(
        adata=adata,
        roi_column="ROI",
        total_cells=config.n_cells,
        random_state=main_config.seed,
        replace=config.replace,
    )
    logger.info("Subsampling complete.")
    logger.info(f"Subsampled dataset contains {len(adata)} cells.")

    # Save anndata object
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")


def subsample_by_roi(
    adata: AnnData,
    roi_column: str = "ROI",
    total_cells: int = 2000,
    random_state: int = 0,
    replace: bool = False,
) -> AnnData:
    """Subsample an AnnData object so that each ROI is represented equally.

    This function selects an equal number of cells from each region of interest (ROI)
    in the input AnnData object. It does not perform any quality control; it only
    performs balanced subsampling.

    Args:
        adata (AnnData):
            Input AnnData object with ROI labels stored in ``adata.obs[roi_column]``.
        roi_column (str):
            Column in ``adata.obs`` indicating ROI labels (e.g., "ROI").
        total_cells (int):
            Desired total number of cells after subsampling. The function calculates
            an equal number per ROI that fits this total.
        random_state (int, optional):
            Seed for reproducible sampling.
        replace (bool, optional):
            If True, sample with replacement (allowing perfect balance even if some
            ROIs are small). If False, the per-ROI sample size is capped by the
            smallest ROIs size.

    Returns:
        AnnData:
            A new AnnData object containing the balanced subsample.

    Notes:
        - The function assumes that each ROI is represented in the dataset.
        - When ``replace=False``, the per-ROI sample size cannot exceed the size of
        the smallest ROI.
    """
    if roi_column not in adata.obs:
        raise KeyError(f"Column '{roi_column}' not found in adata.obs")

    rng = np.random.default_rng(random_state)

    # Count cells per ROI
    roi_counts = adata.obs[roi_column].value_counts()
    rois = roi_counts.index.tolist()
    num_rois = len(rois)
    if num_rois == 0:
        raise ValueError("No ROIs found to stratify by.")

    # Target equal count per ROI based on requested total
    requested_per_roi = total_cells // num_rois
    if requested_per_roi == 0:
        raise ValueError(
            f"total_cells={total_cells} is too small for {num_rois} ROIs. "
            "Increase total_cells so that at least 1 cell per ROI can be sampled."
        )

    # If sampling without replacement, we cannot take more than the smallest ROI size
    if not replace:
        max_per_roi_no_replace = int(roi_counts.min())
        per_roi = min(requested_per_roi, max_per_roi_no_replace)
    else:
        per_roi = requested_per_roi

    # If per_roi is zero after capping, give a helpful error
    if per_roi == 0:
        smallest = int(roi_counts.min())
        raise ValueError(
            "Cannot draw at least 1 cell per ROI without replacement. "
            f"Smallest ROI has only {smallest} cells. Either set replace=True "
            f"or reduce the number of ROIs or adjust total_cells."
        )

    # Collect indices sampled per ROI
    sampled_obs_names = []
    for roi in rois:
        # Get the observation names for this ROI
        roi_mask = (adata.obs[roi_column] == roi).values
        roi_obs_names = adata.obs_names[roi_mask]

        # Draw per_roi samples for this ROI
        # Note: use rng.choice for reproducibility
        chosen = rng.choice(roi_obs_names, size=per_roi, replace=replace)
        sampled_obs_names.extend(chosen.tolist())

    # Subset the AnnData to the sampled cells (preserve order of obs_names)
    sampled_obs_names = pd.Index(sampled_obs_names)
    adata_sub = adata[sampled_obs_names, :].copy()

    # Helpful provenance note
    adata_sub.uns = dict(adata_sub.uns) if adata_sub.uns is not None else {}
    adata_sub.uns["subsample_by_roi"] = {
        "roi_column": roi_column,
        "total_cells_requested": total_cells,
        "num_rois": num_rois,
        "per_roi": per_roi,
        "replace": replace,
        "random_state": random_state,
    }

    return adata_sub
