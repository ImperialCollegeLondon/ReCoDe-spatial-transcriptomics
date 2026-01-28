"""Annotation module."""

import warnings
from logging import getLogger

# import torch
import pandas as pd
import scanpy as sc
import seaborn as sns
import squidpy as sq

from recode_st.config import AnnotateModuleConfig, IOConfig

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)


def remove_small_clusters(adata: sc.AnnData, cluster_name: str) -> sc.AnnData:
    """Remove clusters with less than 10 cells.

    Args:
        adata: AnnData object.
        cluster_name: Name of the cluster column in adata.obs.

    Returns:
        Filtered AnnData object.
    """
    # Remove cell clusters with less than 10 cells
    logger.info("Removing clusters with fewer than 10 cells...")
    cluster_counts = adata.obs[cluster_name].value_counts()
    clusters_to_remove = cluster_counts[cluster_counts < 10].index
    adata = adata[~adata.obs[cluster_name].isin(clusters_to_remove)].copy()
    logger.info(f"Clusters removed after filtering: {clusters_to_remove.tolist()}")
    logger.info(f"Remaining clusters: {adata.obs[cluster_name].unique().tolist()}")


def run_annotate(config: AnnotateModuleConfig, io_config: IOConfig):
    """Run annotation on Xenium data."""
    # Set variables
    module_dir = io_config.output_dir / config.module_name
    cluster_name = config.cluster_name
    new_clusters = config.clusters_label

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set figure directory for this module (overrides global setting)
    sc.settings.figdir = module_dir

    # Set figure settings to ensure consistency across all modules
    cmap = sns.color_palette("Spectral", as_cmap=True)

    # Import data
    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(io_config.output_dir / "dimension_reduction" / "adata.h5ad")

    # Ensure categorical
    adata.obs[cluster_name] = adata.obs[cluster_name].astype("category")
    # Create a palette for the clusters
    color_palette = sns.color_palette(
        "hls", len(adata.obs[cluster_name].cat.categories)
    )
    adata.uns[f"{cluster_name}_colors"] = [color for color in color_palette.as_hex()]
    logger.info(
        f"Saved palette {cluster_name}_colors: {adata.uns[f'{cluster_name}_colors']}"
    )

    # Calculate the differentially expressed genes for every cluster,
    # compared to the rest of the cells in our adata
    logger.info(
        f"Calculating differentially expressed genes for each cluster: {cluster_name}"
    )
    sc.tl.rank_genes_groups(adata, groupby=cluster_name, method="wilcoxon")

    logger.info("Plotting the top differentially expressed genes for each cluster...")
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=cluster_name,
        standard_scale="var",
        n_genes=5,
        cmap=cmap,
        show=False,
        save=f"{config.module_name}_{cluster_name}.pdf",
    )
    logger.info(f"Dotplot saved to {sc.settings.figdir}")

    # Plot differentially expressed genes for each cluster
    logger.info("Plot differentially expressed genes for each cluster in elbow plot...")
    sc.pl.rank_genes_groups(
        adata,
        n_genes=10,
        ncols=3,
        legend_fontsize=10,
        cmap=cmap,
        show=False,
        save=f"_{config.module_name}_{cluster_name}.pdf",
    )
    logger.info(f"UMAP plot saved to {sc.settings.figdir}")

    # Make a dataframe of marker expression
    logger.info("Save files for differentially expressed genes for each cluster...")
    logger.info("File 1...")
    markers = sc.get.rank_genes_groups_df(adata, None)
    markers = markers[(markers["pvals_adj"] < 0.05) & (markers["logfoldchanges"] > 0.5)]
    markers.to_excel(
        module_dir / "markers.xlsx",
        index=False,
    )
    logger.info(f"Markers saved to {module_dir}")

    logger.info("File 2...")
    # Define the number of clusters
    clusters_list = len(adata.obs[cluster_name].astype(str).unique())

    # Create a list
    list = []
    for cluster_number in range(clusters_list):
        top_genes = adata.uns["rank_genes_groups"]["names"][
            str(cluster_number)
        ]  # Get the names of the top differentially expressed genes
        top_genes = top_genes[:10]  # Get the top 10 genes
        new_row = pd.Series(
            {"Cluster Number": cluster_number, "Top Genes": top_genes}
        )  # Create a new row with the cluster_number and top_genes
        list.append(new_row)

    # Convert list of series to DataFrame
    diff_gene_df = pd.concat(list, axis=1).T
    diff_gene_df.set_index(diff_gene_df.columns[0], inplace=True)
    diff_gene_df.to_csv(
        module_dir / "top_differentially_expressed_genes.csv",
        index=True,
    )
    logger.info(f"Top differentially expressed genes saved to {module_dir}")

    logger.info("File 3...")
    # Create a dictionary to store DataFrames for each cluster
    cluster_dict = {}
    cluster_path = module_dir / "cluster_diff_genes"
    cluster_path.mkdir(exist_ok=True)
    for cluster_number in range(clusters_list):
        # print(cluster_number)
        current_cluster = markers[markers["group"] == str(cluster_number)].sort_values(
            by="logfoldchanges", ascending=False
        )  # make a dataframe of the current cluster
        cluster_dict[f"cluster_{cluster_number}"] = (
            current_cluster  # Store the DataFrame in the dictionary
        )
        # Export the DataFrame to a CSV file
        csv_filename = cluster_path / f"cluster_{cluster_number}_data.csv"

        current_cluster.to_csv(csv_filename, index=False)
        logger.info(f"Exported cluster {cluster_number} data to {csv_filename}")

    # Rename the clusters based on the markers
    cluster_to_cell_type = config.cluster_to_cell_type
    logger.info(f"Cluster to cell type mapping: {cluster_to_cell_type}")
    logger.info(f"Type of data in mapping: {type(cluster_to_cell_type)}")

    logger.info("Renaming clusters based on markers...")
    cluster_to_cell_type_dict = config.cluster_to_cell_type  # import from config
    adata.obs[new_clusters] = adata.obs[cluster_name].map(cluster_to_cell_type_dict)

    # Ensure categorical
    adata.obs[new_clusters] = adata.obs[new_clusters].astype("category")
    # Create a palette for the new clusters
    color_palette = sns.color_palette("hls", len(adata.obs[new_clusters].unique()))
    adata.uns[f"{new_clusters}_colors"] = [color for color in color_palette.as_hex()]
    logger.info(
        f"Saved palette {new_clusters}_colors: {adata.uns[f'{new_clusters}_colors']}"
    )

    logger.info("Plotting UMAP with new cluster names...")
    sc.pl.umap(
        adata,
        color=[cluster_name, new_clusters],
        legend_loc="right margin",
        legend_fontsize=12,
        frameon=False,
        ncols=2,  # Side by side
        wspace=0.4,  # Space between plots
        title=new_clusters,
        show=False,
        save=f"_{config.module_name}_{new_clusters}_combined_annotation.pdf",
    )
    logger.info(f"UMAP plot with new cluster names saved to {sc.settings.figdir}")

    # Calculate the differentially expressed genes for every cluster,
    # compared to the rest of the cells in our adata
    logger.info(
        f"Calculating differentially expressed genes for each cluster: {new_clusters}"
    )
    sc.tl.rank_genes_groups(adata, groupby=new_clusters, method="wilcoxon")

    logger.info("Plotting the top differentially expressed genes for each cluster...")
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=new_clusters,
        standard_scale="var",
        n_genes=5,
        cmap=cmap,
        show=False,
        save=f"{config.module_name}_{new_clusters}.pdf",
    )
    logger.info(f"Dotplot saved to {sc.settings.figdir}")

    # View specific gene expression
    logger.info("Plotting genes of interest on tissue...")
    ROI_list = adata.obs["ROI"].unique().tolist()
    for roi in ROI_list:
        adata_roi = adata[adata.obs["ROI"] == roi]
        sq.pl.spatial_scatter(
            adata_roi,
            library_id="spatial",
            shape=None,
            outline=False,
            color=new_clusters,
            size=0.5,
            figsize=(15, 15),
            save=f"clusters_{roi}.pdf",
            dpi=300,
        )
        logger.info(f"Saved cluster plot for ROI {roi} to {module_dir}")

    # Save anndata object
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")
    logger.info("Annotation module completed successfully.")
