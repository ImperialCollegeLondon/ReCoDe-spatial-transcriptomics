"""Calculate pseudobulk differential expression for Xenium data."""

import warnings
from logging import getLogger

import matplotlib.pyplot as plt
import numpy as np

# import torch
import pandas as pd
import scanpy as sc
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from recode_st.config import IOConfig, PsuedobulkModuleConfig

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)


# Functions
def filter_small_rois(adata: sc.AnnData, min_cells: int = 100):
    """Remove ROIs with fewer than min_cells cells.

    Args:
        adata: AnnData
        min_cells: Minimum number of cells per ROI.

    Returns:
        Filtered AnnData.
    """
    roi_counts = adata.obs["ROI"].value_counts()
    rois_to_remove = roi_counts[roi_counts < min_cells].index

    if len(rois_to_remove) > 0:
        logger.info("Removed ROIs:", list(rois_to_remove))
    else:
        logger.info("All ROIs have at least", min_cells, "cells. Nothing removed.")

    return adata[~adata.obs["ROI"].isin(rois_to_remove)].copy()


def volcano_plot(
    de_results: pd.DataFrame,
    lfc_threshold: int = 5,
    padj_threshold: float = 0.05,
    use_padj=True,
    title="Volcano Plot: Differential Expression",
):
    """Volcano plot for differentially expressed genes.

    Args:
        de_results : DeseqStats or pd.DataFrame DESeq2 results object or DataFrame.
        lfc_threshold : float
        Absolute log2 fold change threshold to highlight genes.
        padj_threshold : float
        Adjusted p-value threshold to highlight genes.
        use_padj : bool
        If True, use adjusted p-value (padj)
        title : str
        Plot title.

    Returns:
        Plot.
    """
    # Convert to DataFrame if DeseqStats
    if hasattr(de_results, "results_df"):
        df = de_results.results_df.copy()
    else:
        df = de_results.copy()

    # Determine p-value column
    p_col = "padj"
    if p_col not in df.columns:
        raise ValueError(f"{p_col} column not found in DE results")

    # Fill NaNs
    df[p_col] = df[p_col].fillna(1)

    # Determine gene categories
    up = (df["log2FoldChange"] > lfc_threshold) & (df[p_col] < padj_threshold)
    down = (df["log2FoldChange"] < -lfc_threshold) & (df[p_col] < padj_threshold)
    sig_small_lfc = (abs(df["log2FoldChange"]) <= lfc_threshold) & (
        df[p_col] < padj_threshold
    )
    neutral = ~(up | down | sig_small_lfc)

    plt.figure(figsize=(10, 12))

    # Plot categories
    plt.scatter(
        df.loc[neutral, "log2FoldChange"],
        -np.log10(df.loc[neutral, p_col]),
        color="grey",
        alpha=0.6,
        s=10,
        fontsize=8,
    )
    plt.scatter(
        df.loc[up, "log2FoldChange"],
        -np.log10(df.loc[up, p_col]),
        color="#8ab184",
        alpha=0.8,
        label=f"Upregulated (LFC > {lfc_threshold})",
        s=15,
    )
    plt.scatter(
        df.loc[down, "log2FoldChange"],
        -np.log10(df.loc[down, p_col]),
        color="#BC3C29",
        alpha=0.8,
        label=f"Downregulated (LFC < -{lfc_threshold})",
        s=15,
    )
    plt.scatter(
        df.loc[sig_small_lfc, "log2FoldChange"],
        -np.log10(df.loc[sig_small_lfc, p_col]),
        color="#6F99AD",
        alpha=0.8,
        label=f"Significant, |LFC| ≤ {lfc_threshold}",
        s=15,
    )

    # Label only up/down genes
    for _, row in df.loc[up | down].iterrows():
        plt.text(
            row["log2FoldChange"],
            -np.log10(row[p_col]),
            row.name,
            fontsize=12,
            ha="right",
        )

    # Add vertical dotted lines for LFC threshold
    plt.axvline(lfc_threshold, color="black", linestyle=":", linewidth=1)
    plt.axvline(-lfc_threshold, color="black", linestyle=":", linewidth=1)

    # Add horizontal dotted line for p-value threshold
    plt.axhline(-np.log10(padj_threshold), color="black", linestyle=":", linewidth=1)

    # Labels
    plt.xlabel("log2 Fold Change (COPD vs IPF)")
    plt.ylabel(f"-log10({p_col})")
    plt.title(title)
    plt.legend()
    plt.tight_layout()


def run_annotate(config: PsuedobulkModuleConfig, io_config: IOConfig):
    """Run annotation on Xenium data."""
    # Set variables
    module_dir = io_config.output_dir / config.module_name
    # Define parameters
    res = 0.5
    subset = "immune"
    annotation_col = subset + "_annotation"
    subset_key = f"{subset}_{res}"

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set figure directory for this module (overrides global setting)
    sc.settings.figdir = module_dir

    # Set figure settings to ensure consistency across all modules
    cmap = sns.color_palette("Spectral", as_cmap=True)

    # Import data
    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(io_config.output_dir / "dimension_reduction" / "adata.h5ad")
    adata = sc.read_h5ad(f"manual/subset/{subset}_adata.h5ad")
    logger.info("Xenium data loaded.")

    # Subset to only desired annotation
    logger.info(
        f"Calculating DEG for each cell type in {annotation_col} annotations..."
    )
    unique_cell_types = adata.obs[annotation_col].unique().tolist()
    for cell_type in unique_cell_types:
        logger.info("Processing cell type:", cell_type)
        cell_subset = adata[adata.obs[annotation_col] == cell_type]

        # Ensure all samples have at least 100 cells
        cell_subset = filter_small_rois(cell_subset, min_cells=50)

        pbs = []  # pseudobulk sample
        for sample in cell_subset.obs["ROI"].unique():
            sample_cell_subset = cell_subset[cell_subset.obs["ROI"] == sample].copy()
            sample_cell_subset.X = sample_cell_subset.layers[
                "counts"
            ]  # Use raw data rather than normalized log transformed data
            # Create new adata for each sample
            rep_data = sc.AnnData(
                X=sample_cell_subset.X.sum(axis=0), var=sample_cell_subset.var[[]]
            )

            rep_data.obs_names = [sample]  # set obs name to the sample
            rep_data.obs["condition"] = sample_cell_subset.obs["condition"].iloc[0]

            # Add to psuedobulk
            pbs.append(rep_data)

            # Concatenate all pseudobulk samples
            psuedobulk = sc.concat(pbs)

            # Differential expression
            counts = pd.DataFrame(
                psuedobulk.X, columns=psuedobulk.var_names
            )  # need to do this to pass var names

            # Create DESeq2 dataset
            dds = DeseqDataSet(
                counts=counts, metadata=psuedobulk.obs, design_factors="condition"
            )

            # Filter to remove genes not found in at least one sample
            sc.pp.filter_genes(dds, min_cells=1)
            print(f"Number of genes after filtering: {dds.X.shape[1]}")

            # Initialize DESeq2
            dds.deseq2()

            # Get DE results
            de_results = DeseqStats(
                dds, n_cpus=8, contrast=("condition", "COPD", "IPF")
            )
            de_results.summary()
            results_df = de_results.results_df  # Get results dataframe
            results_df.to_csv(
                f"pseudobulk_DE_results_{cell_type}.csv"
            )  # Save DE results

            # Plot PCA of psuedobulk samples
            sc.tl.pca(dds)
            sc.pl.pca(
                dds,
                color="condition",
                size=200,
                title=f"PCA of psuedobulk samples - {cell_type}",
                save=f"_pseudobulk_PCA_{cell_type}.pdf",
            )

            # Plot volcano plot
            volcano_plot(
                de_results, lfc_threshold=2, padj_threshold=0.05, use_padj=True
            )
            plt.savefig(f"volcano_plot_pseudobulk_{cell_type}.pdf")
            plt.close()
    logger.info("Pseudobulk differential expression complete.")
