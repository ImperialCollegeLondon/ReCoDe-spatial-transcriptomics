"""Drug2Cell module."""

import warnings
from logging import getLogger

import drug2cell as d2c
import scanpy as sc
import seaborn as sns
import squidpy as sq

from recode_st.config import Drug2CellModuleConfig, IOConfig
from recode_st.helper_function import configure_scanpy_figures

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)

# Variables
CELL_TYPE_TOP = "ingest_pred_cell_type"
CELL_TYPE_LEVEL = "transf_ann_level_2_label"
KEY_CELL_TYPE = "Lymphoid"


def drug2cell_calculation(adata, module_dir):
    """Calculate drug2cell scores.

    Args:
        adata (anndata.AnnData): Annotated data object
        module_dir (Path): Directory to save module outputs

    Returns:
        anndata.AnnData or None: adata with drug2cell scores added, or None if failed
    """
    # Computes the mean of the expression of each gene group in each cell.
    # By default, the function will load a set of ChEMBL drugs and their targets
    drug2cell_path = module_dir / "adata_drug2cell.h5ad"

    if drug2cell_path.exists():
        logger.info("Loading existing drug2cell results...")
        adata = sc.read_h5ad(module_dir / "adata_drug2cell.h5ad")
        return adata
    else:
        logger.info("Calculating drug2cell scores...")
        try:
            d2c.score(adata, use_raw=False)
        except Exception as e:
            logger.error(f"Failed to calculate drug2cell scores: {e}")
            return None

    # Validate that drug2cell results were generated
    if "drug2cell" not in adata.uns:
        logger.error(
            "drug2cell scores not found in adata.uns. Drug scoring may have failed."
        )
        return None

    logger.info(
        f"Successfully calculated drug scores for "
        f"{adata.uns['drug2cell'].shape[1]} drugs"
    )

    # Save results
    try:
        logger.info("Saving drug list CSV...")
        drugs_present = adata.uns["drug2cell"].var
        drugs_present.to_csv(module_dir / "drug2cell_drugs.csv")

        logger.info("Saving adata with drug2cell scores...")
        adata.write_h5ad(drug2cell_path)
        logger.info(f"Results saved to {drug2cell_path}")
    except Exception as e:
        logger.warning(f"Failed to save results: {e}. Continuing anyway...")

    return adata


def visualize_drug2cell(config, adata, cell_type_top: str = CELL_TYPE_TOP, cmap=None):
    """Visualize drug2cell scores.

    Args:
        config (Drug2CellModuleConfig): Configuration object
        adata (anndata.AnnData): Annotated data object with drug2cell scores
        cell_type_top (str): Column name for cell type annotations
        cmap: Color map for visualization

    Returns:
        None
    """
    logger.info("Checking drugs in drug_list against calculated drugs...")
    drug_list = config.drug_list
    available_drugs = adata.uns["drug2cell"].var_names.tolist()
    missing_drugs = [drug for drug in drug_list if drug not in available_drugs]
    if missing_drugs:
        logger.warning(
            f"The following drugs from the config are not found in the "
            f"calculated drug2cell scores: {missing_drugs}"
        )
        drug_list = [drug for drug in drug_list if drug in available_drugs]
        if not drug_list:
            logger.error("No valid drugs found for visualization. Skipping.")
            return

    logger.info("Visualize drug score in UMAP")
    sc.pl.umap(
        adata.uns["drug2cell"],
        color=drug_list,
        color_map=cmap,
        save=f"_{config.module_name}_drug2cell_umap.png",
    )
    logger.info("Visualize drug score in spatial scatter plots")
    for roi in adata.obs["ROI"].unique():
        logger.info(f"Visualizing ROI: {roi}")
        subset = adata[adata.obs["ROI"] == roi]
        # Copy and subset drug2cell scores to match the ROI subset
        subset.uns["drug2cell"] = adata.uns["drug2cell"][subset.obs_names, :]
        sq.pl.spatial_scatter(
            subset.uns["drug2cell"],
            library_id="spatial",
            size=0.5,
            shape=None,
            color=drug_list,
            wspace=0.4,
            cmap=cmap,
            save=f"_{config.module_name}_{roi}_scatter.png",
        )


def filter_low_count_groups(
    adata, cell_type_top: str = CELL_TYPE_TOP, min_cells_per_group: int = 10
):
    """Filter adata to only include groups with sufficient cell counts.

    Args:
        adata (anndata.AnnData): Annotated data object
        cell_type_top (str, optional):Defaults to CELL_TYPE_TOP.
        min_cells_per_group (int, optional): Minimum number of cells per group.

    Returns:
        adata_filtered (): filtered adata object
    """
    # Check if we have enough cells per group for differential expression
    group_counts = adata.obs[cell_type_top].value_counts()
    valid_groups = group_counts[group_counts >= min_cells_per_group].index.tolist()

    if len(valid_groups) < 2:
        logger.warning(
            f"Not enough groups with >= {min_cells_per_group} cells "
            f"Skipping filtering adata object."
        )
    else:
        # Filter to only valid groups
        valid_mask = adata.obs[cell_type_top].isin(valid_groups)
        adata_filtered = adata.uns["drug2cell"][valid_mask, :]

        logger.info(
            f"Filtered for {len(valid_groups)} groups "
            f"with >= {min_cells_per_group} cells each"
        )
        return adata_filtered


def calc_DE(adata, config, cell_type_top: str = CELL_TYPE_TOP, cmap=None):
    """Calculate differential expression for drug2cell scores.

    Args:
        adata (anndata.AnnData): Annotated data object
        config (Drug2CellModuleConfig): Configuration object
        cell_type_top (str): Column name for cell type annotations
        cmap: Color map for visualization
    Returns:
        None
    """
    logger.info("Calculating differential expression for drug2cell scores...")
    sc.tl.rank_genes_groups(adata, method="wilcoxon", groupby=CELL_TYPE_TOP)

    logger.info("Visualizing differential expression results...")
    sc.pl.rank_genes_groups_dotplot(
        adata,
        swap_axes=True,
        dendrogram=False,
        n_genes=3,
        cmap=cmap,
        save=f"_{config.module_name}_{CELL_TYPE_TOP}_dotplot.png",
    )


def celltype_level(
    config,
    adata,
    cell_type: str = KEY_CELL_TYPE,
    cell_type_level: str = CELL_TYPE_LEVEL,
    cell_type_top: str = CELL_TYPE_TOP,
    cmap=None,
):
    """Calculate drug2cell for specified cell type for set cell type level.

    Args:
        config: Configuration object
        adata (anndata.AnnData): Annotated data object
        cell_type (str): Cell type to subset on
        cell_type_level (str): Cell type level or resolution of interest
        cell_type_top (str): Column name for cell type annotations
        cmap: Color map for visualization

    Returns:
        None
    """
    # Fix parameter name mismatch: CLUSTER_TYPE should be CELL_TYPE_LEVEL_ALL
    logger.info("Checking drugs in drug_list against calculated drugs...")
    drug_list = config.drug_list

    logger.info(f"Subsetting on {cell_type} in {cell_type_level}")
    # Fix critical bug: was using assignment (=) instead of comparison (==)
    if cell_type_level not in adata.obs.columns:
        logger.error(f"Column '{cell_type_level}' not found in adata.obs")
        return

    if cell_type not in adata.obs[cell_type_level].values:
        logger.warning(
            f"Cell type '{cell_type}' not found in column '{cell_type_level}'"
        )
        return

    subset = adata[adata.obs[cell_type_level] == cell_type]

    if len(subset) == 0:
        logger.warning(
            f"No cells found for {cell_type} in {cell_type_level}. "
            "Skipping visualization."
        )
        return

    # Copy drug2cell scores to subset if they exist
    if "drug2cell" in adata.uns:
        # Properly subset drug2cell scores to match only the myeloid cells
        subset.uns["drug2cell"] = adata.uns["drug2cell"][subset.obs_names, :]
    else:
        logger.error("No drug2cell scores found in adata.uns. Cannot visualize subset.")
        return

    logger.info(f"Visualizing on {cell_type} in {cell_type_level}")

    logger.info("Visualize drug score in UMAP")
    sc.pl.umap(
        subset.uns["drug2cell"],
        color=[*drug_list, cell_type_top],
        color_map=cmap,
        save=f"_{config.module_name}_drug2cell_{cell_type}.png",
    )

    logger.info("Visualize drug score with dotplot")
    sc.pl.dotplot(
        subset.uns["drug2cell"],
        swap_axes=True,
        cmap=cmap,
        n_genes=5,
        save=f"_{config.module_name}_{cell_type_top}_{cell_type}.png",
    )


def run_drug2cell(config: Drug2CellModuleConfig, io_config: IOConfig):
    """Calculate drug score using drug2cell.

    Args:
        config (Drug2CellModuleConfig): Drug2Cell module configuration
        io_config (IOConfig): IO configuration object

    Returns:
        None
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

    logger.info("Starting Drug2Cell module...")

    logger.info("Loading Xenium data...")
    # Fix path - should load from integrate_ingest module output
    input_file = io_config.output_dir / "3_integrate_ingest" / "adata.h5ad"
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        logger.error("Please run the integrate_ingest module first.")
        return

    adata = sc.read_h5ad(input_file)

    logger.info("Calculating drug score every cell in adata using ChEMBL database")
    adata = drug2cell_calculation(adata, module_dir)
    if adata is None:
        logger.error("Failed to calculate drug2cell scores. Aborting module.")
        return

    logger.info("Filtering low count groups for differential expression...")
    adata = filter_low_count_groups(adata, cell_type_top=CELL_TYPE_TOP)
    if adata is None:
        logger.error("Failed to filter low count groups. Aborting module.")
        return

    logger.info("Visualize drug2cell...")
    # Validate required columns exist
    if CELL_TYPE_LEVEL not in adata.obs.columns:
        logger.warning(
            f"Column '{CELL_TYPE_LEVEL}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )
        logger.info("Proceeding with cell_type_level_all only...")
    visualize_drug2cell(config, adata=adata, cell_type_top=CELL_TYPE_TOP, cmap=cmap)

    logger.info("Calculating differential expression...")
    calc_DE(adata, config, cell_type_top=CELL_TYPE_TOP, cmap=cmap)

    logger.info("Visualize only respiratory drugs")
    # Apply the same filtering as used in differential expression
    group_counts = adata.obs[CELL_TYPE_TOP].value_counts()
    min_cells_per_group = 10  # Same threshold as calc_DE
    valid_groups = group_counts[group_counts >= min_cells_per_group].index.tolist()

    if len(valid_groups) >= 2:
        # Filter to only valid groups for respiratory drugs visualization
        drug2cell_adata = adata.uns["drug2cell"].copy()
        valid_mask = drug2cell_adata.obs[CELL_TYPE_TOP].isin(valid_groups)
        drug2cell_filtered = drug2cell_adata[valid_mask, :].copy()

        plot_args = d2c.util.prepare_plot_args(drug2cell_filtered, categories=["R"])
        sc.pl.dotplot(
            drug2cell_filtered,
            groupby=CELL_TYPE_TOP,
            swap_axes=True,
            **plot_args,
            cmap=cmap,
            save=f"_{config.module_name}_{CELL_TYPE_TOP}_respiratory_drugs_by.png",
        )
    else:
        logger.warning(
            f"Not enough groups with >= {min_cells_per_group} cells "
            "for respiratory drugs visualization. Skipping."
        )

    # Only run cell type specific analysis if the column exists
    if CELL_TYPE_LEVEL in adata.obs.columns:
        logger.info(f"Examine drug score in {KEY_CELL_TYPE} in {CELL_TYPE_LEVEL}")
        celltype_level(
            config, adata, KEY_CELL_TYPE, CELL_TYPE_LEVEL, CELL_TYPE_TOP, cmap
        )
    else:
        logger.info(
            f"Skipping cell type specific analysis - "
            f"column '{CELL_TYPE_LEVEL}' not available"
        )

    logger.info("Drug2Cell module finished.")
