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
ANNOTATION = "manual_annotation"


def drug2cell_calculation(adata, module_dir):
    """Calculate drug2cell scores.

    Args:
        adata (anndata.AnnData): Annotated data object
        module_dir (Path): Directory to save module outputs

    Returns:
        AnnData or None: adata with drug2cell scores added, or None if failed
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


def filter_low_count_groups(
    adata, annotation: str = ANNOTATION, min_cells_per_group: int = 10
):
    """Filter adata to only include groups with sufficient cell counts.

    Args:
        adata (anndata.AnnData): Annotated data object
        annotation (str, optional):Defaults to ANNOTATION.
        min_cells_per_group (int, optional): Minimum number of cells per group.

    Returns:
        adata_filtered (): filtered adata object
    """
    # Check if we have enough cells per group for differential expression
    group_counts = adata.obs[annotation].value_counts()
    valid_groups = group_counts[group_counts >= min_cells_per_group].index.tolist()

    if len(valid_groups) < 2:
        logger.warning(
            f"Not enough groups with >= {min_cells_per_group} cells "
            f"Skipping filtering adata object."
        )
    else:
        valid_mask = adata.obs[annotation].isin(valid_groups)
        adata_filtered = adata.uns["drug2cell"][valid_mask, :]

        logger.info(
            f"Filtered for {len(valid_groups)} groups "
            f"with >= {min_cells_per_group} cells each"
        )
        return adata_filtered


def visualize_drug2cell(config, adata, annotation: str = ANNOTATION, cmap=None):
    """Visualize drug2cell scores.

    Args:
        config (Drug2CellModuleConfig): Configuration object
        adata (anndata.AnnData): Annotated data object with drug2cell scores
        annotation (str): Column name for cell type annotations
        cmap: Color map for visualization

    Returns:
        None
    """
    logger.info("Checking drugs in drug_list against calculated drugs...")

    # Check if drug2cell results exist
    if "drug2cell" not in adata.uns:
        logger.error(
            "drug2cell scores not found in adata.uns. "
            "Drug scoring may have failed or not been performed."
        )
        logger.info("Skipping drug2cell visualization.")
        return

    drug_list = config.drug_list  # the drugs you want to visualize
    available_drugs = adata.uns[
        "drug2cell"
    ].var_names.tolist()  # all drugs with calculated scores
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
        save=f"_{config.module_name}_drug2cell_umap.pdf",
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
            save=f"_{config.module_name}_{roi}_scatter.pdf",
        )


def calc_n_plot_DE(adata, config, annotation: str = ANNOTATION, cmap=None):
    """Calculate differential expression for drug2cell scores.

    Args:
        adata (anndata.AnnData): Annotated data object
        config (Drug2CellModuleConfig): Configuration object
        annotation (str): Column name for cell type annotations
        cmap: Color map for visualization
    Returns:
        None
    """
    logger.info("Calculating differential expression for drug2cell scores...")
    sc.tl.rank_genes_groups(
        adata.uns["drug2cell"], method="wilcoxon", groupby=ANNOTATION
    )

    logger.info("Visualizing differential expression results...")
    sc.pl.rank_genes_groups_dotplot(
        adata.uns["drug2cell"],
        swap_axes=True,
        dendrogram=False,
        n_genes=3,
        cmap=cmap,
        save=f"_{config.module_name}_{ANNOTATION}_dotplot.pdf",
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
    cmap_blues = sns.color_palette("Blues", as_cmap=True)

    logger.info("Starting Drug2Cell module...")

    logger.info("Loading Xenium data...")
    input_file = io_config.output_dir / "annotate" / "adata.h5ad"
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        logger.error("Please ensure correct file path.")
        return

    adata = sc.read_h5ad(input_file)

    logger.info("Calculating drug score every cell in adata using ChEMBL database")
    adata = drug2cell_calculation(adata, module_dir)
    if adata is None:
        logger.error("Failed to calculate drug2cell scores. Aborting module.")
        return

    logger.info("Visualize dotplot of drug2cell scores across cell types")
    sc.pl.dotplot(
        adata.uns["drug2cell"],
        var_names=config.drug_list,
        groupby=ANNOTATION,
        swap_axes=True,
        cmap=cmap,
        save=f"_{config.module_name}_{ANNOTATION}_all_drugs_by.pdf",
    )

    logger.info("Filtering low count groups for differential expression...")
    filtered_drug2cell = filter_low_count_groups(adata, annotation=ANNOTATION)
    if filtered_drug2cell is None:
        logger.error("Failed to filter low count groups. Aborting module.")
        return

    logger.info("Visualize drug2cell...")
    # Validate required columns exist
    if ANNOTATION not in adata.obs.columns:
        logger.warning(
            f"Column '{ANNOTATION}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )
    else:
        # Visualize drugs from drug list in a UMAP (all) and spatial scatter (each ROI)
        visualize_drug2cell(config, adata=adata, annotation=ANNOTATION, cmap=cmap_blues)

    if "drug2cell" in adata.uns:
        logger.info("Calculating differential expression...")
        calc_n_plot_DE(adata, config, annotation=ANNOTATION, cmap=cmap)
    else:
        logger.warning(
            "Skipping differential expression analysis - no drug2cell results available"
        )

    logger.info("Visualize only respiratory drugs")

    # Check if drug2cell results are available
    if "drug2cell" not in adata.uns:
        logger.warning(
            "Skipping respiratory drugs visualization - no drug2cell results available"
        )
        logger.info("Drug2Cell module finished.")
        return

    # Apply the same filtering as used in differential expression
    group_counts = adata.obs[ANNOTATION].value_counts()
    min_cells_per_group = 10  # Same threshold as calc_DE
    valid_groups = group_counts[group_counts >= min_cells_per_group].index.tolist()

    if len(valid_groups) >= 2:
        # Filter to only valid groups for respiratory drugs visualization
        drug2cell_adata = adata.uns["drug2cell"].copy()
        valid_mask = drug2cell_adata.obs[ANNOTATION].isin(valid_groups)
        drug2cell_filtered = drug2cell_adata[valid_mask, :].copy()

        plot_args = d2c.util.prepare_plot_args(drug2cell_filtered, categories=["R"])
        sc.pl.dotplot(
            drug2cell_filtered,
            groupby=ANNOTATION,
            swap_axes=True,
            **plot_args,
            cmap=cmap,
            save=f"_{config.module_name}_{ANNOTATION}_respiratory_drugs_by.pdf",
        )
    else:
        logger.warning(
            f"Not enough groups with >= {min_cells_per_group} cells "
            "for respiratory drugs visualization. Skipping."
        )

    logger.info("Drug2Cell module finished.")
