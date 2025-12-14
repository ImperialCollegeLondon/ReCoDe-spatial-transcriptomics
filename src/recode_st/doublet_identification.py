"""Identify doublets and spatial overlap with ovrly."""

import pickle
import warnings
from logging import getLogger
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ovrlpy
import polars as pl
import scanpy as sc
import seaborn as sns

from recode_st.config import DoubletIdentificationModuleConfig, IOConfig
from recode_st.helper_function import configure_scanpy_figures

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


def load_transcripts(io_config: IOConfig):
    """Load transcripts file for each ROI."""
    xenium_path = Path(io_config.xenium_dir)
    output_data_dir = Path(io_config.output_data_dir)
    output_data_dir.mkdir(parents=True, exist_ok=True)

    if not xenium_path.exists():
        raise FileNotFoundError(f"Xenium root directory not found: {xenium_path}")

    all_dfs = {}
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
                    logger.info(f"Loading Xenium transcript data for {roi_name}...")
                    df = ovrlpy.io.read_Xenium(roi_folder / "transcripts.parquet")
                    all_dfs[roi_name] = df
                    logger.info(
                        f"Loaded Xenium transcript data for {roi_name} saved to dict."
                    )
                except Exception as err:
                    logger.error(f"Failed loading Xenium data for {roi_name}: {err}")
                    continue

    logger.info(f"Loaded transcripts for {len(all_dfs)} ROIs: {list(all_dfs.keys())}")
    return all_dfs


def process_roi_doublets(roi_name, roi_df, module_dir, cmap):
    """Process doublet analysis for a single ROI.

    Args:
        roi_name (str): Name of the ROI
        roi_df (polars.DataFrame): Transcript data for the ROI
        module_dir (Path): Output directory for saving results
        cmap: Colormap for visualization
    """
    logger.info(f"Processing doublets for ROI: {roi_name}")
    logger.info(f"Creating new directory for {roi_name} results...")
    roi_output_dir = module_dir / roi_name
    roi_output_dir.mkdir(exist_ok=True)

    logger.info(f"Identifying doublets for {roi_name} using ovrly ...")
    n = 100
    total_transcripts = roi_df.shape[0]
    transcripts_shown = total_transcripts // n
    percent_shown = (transcripts_shown / total_transcripts) * 100
    logger.info(f"Showing every {n}th transcript for {roi_name}: ")
    logger.info(
        f"Total transcripts: {total_transcripts}, "
        f"Shown: {transcripts_shown} "
        f"({percent_shown:.2f}%)"
    )

    # Plot transcripts
    fig, ax = plt.subplots()
    ax.scatter(roi_df[::n, "x"], roi_df[::n, "y"], s=0.1)
    _ = ax.set(aspect="equal")
    plt.savefig(
        roi_output_dir / f"{roi_name}_transcripts_overview.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()

    # TODO: Subset on HVG? Cell type specific genes? NEED TO THINK

    logger.info("Running the ovrlpy pipeline...")
    roi_ob = ovrlpy.Ovrlp(
        roi_df,
        n_workers=16,  # I think this is the number of CPUs used?
        random_state=42,  # TODO: Determine if I need this since seed
        # everything should deal with this
    )
    roi_ob.analyse()

    logger.info("Visualizing results ...")
    fig = ovrlpy.plot_pseudocells(roi_ob)
    plt.savefig(
        roi_output_dir / f"{roi_name}_pseudocells.png", dpi=300, bbox_inches="tight"
    )
    plt.close()

    # plot SVI per ROI
    fig = ovrlpy.plot_signal_integrity(roi_ob, signal_threshold=3)
    plt.savefig(
        roi_output_dir / f"{roi_name}_signal_integrity.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()

    logger.info("Detecting doublets...")
    doublets = roi_ob.detect_doublets(min_signal=3, integrity_sigma=2)
    logger.info(f"Detected {len(doublets)} doublets in {roi_name}")

    logger.info("Saving doublets to parquet file...")
    doublets.write_parquet(roi_output_dir / "data.parquet")

    logger.info("Visualizing doublets...")
    # Overview of tissue
    fig, ax = plt.subplots()
    _scatter = ax.scatter(
        doublets["x"], doublets["y"], c=doublets["integrity"], s=0.2, cmap=cmap
    )
    _ = ax.set_aspect("equal")
    _ = fig.colorbar(_scatter, ax=ax)
    plt.savefig(
        roi_output_dir / f"{roi_name}_doublets_overview.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()

    # Visualize specific doublet case
    # select random doublet index
    doublet_case = np.random.randint(0, len(doublets))
    x, y = doublets["x"][doublet_case], doublets["y"][doublet_case]
    _ = ovrlpy.plot_region_of_interest(roi_ob, x, y, window_size=60)
    plt.savefig(
        roi_output_dir / f"{roi_name}_doublets_{doublet_case}.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()

    logger.info("Visualize cells in Z stacks using pickle...")
    # Make folder to save pickle files
    pickle_dir = roi_output_dir / "pickle_doublets"
    pickle_dir.mkdir(exist_ok=True)

    # Name pickle file
    pickle_file = pickle_dir / f"{roi_name}.pickle"

    # Save the analyzed roi_ob for future use
    logger.info(f"Saving analysis to {pickle_file}...")
    with open(pickle_file, "wb") as f:
        pickle.dump(roi_ob, f)

    logger.info("Plotting cells in Z stacks...")
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    plt.title(f"3D visualization of {roi_name}")
    for i in range(-2, 3):
        subset = roi_ob.transcripts.filter(
            (pl.col("z") - pl.col("z_center")).is_between(i, i + 1)
        )
        # downsample the number of transcripts
        subset = subset[::100]
        ax.scatter(subset["x"], subset["y"], i, s=1, alpha=0.1)

    ratio = roi_ob.transcripts["x"].max() / roi_ob.transcripts["y"].max()
    ax.set_box_aspect([ratio, 1, 0.75])
    plt.savefig(
        roi_output_dir / f"{roi_name}_3d_transcripts.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


def run_doublet_id(io_config: IOConfig, config: DoubletIdentificationModuleConfig):
    """Identify doublets using ovrl.

    Args:
        io_config (_type_): _description_
        config (_type_): _description_

    Returns:
        _type_: _description_
    """
    # Name of the column to store label transfer results in adata.obs
    module_dir = io_config.output_dir / config.module_name

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set figure directory for this module (overrides global setting)
    sc.settings.figdir = module_dir

    # Set figure settings to ensure consistency across all modules
    configure_scanpy_figures(str(io_config.output_dir))
    cmap = sns.color_palette("Spectral", as_cmap=True)

    logger.info("Loading transcripts form each ROI into a df...")
    # Load "transcripts.parquet" for each ROI and make a dict of dfs
    df_all = load_transcripts(io_config)

    logger.info("Identifying doublets from each ROI...")
    for roi_name, roi_df in df_all.items():
        process_roi_doublets(roi_name, roi_df, module_dir, cmap)

    logger.info("Completed doublet identification for all ROIs.")
