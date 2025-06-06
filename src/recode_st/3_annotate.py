"""Annotation module."""

import os
import warnings  # ? what is the best way to suppress warnings from package inputs?
from logging import getLogger

# import torch
import pandas as pd
import scanpy as sc

from recode_st.helper_function import seed_everything
from recode_st.logging_config import configure_logging
from recode_st.paths import output_path

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def run_annotate():
    """Run annotation on Xenium data."""
    # Set variables
    module_name = "3_annotate"
    module_dir = output_path / module_name

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Import data
    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(output_path / "2_DR/adata.h5ad")

    # change directory to output_path/module_name
    os.chdir(module_dir)
    logger.info(f"Changed directory to {module_dir}")

    # Annotate cell clusters
    # Calculate the differentially expressed genes for every cluster,
    # compared to the rest of the cells in our adata

    logger.info("Calculating differentially expressed genes for each cluster...")
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

    logger.info("Plotting the top differentially expressed genes for each cluster...")
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby="leiden",
        standard_scale="var",
        n_genes=5,
        show=False,
        save="dotplot.png",
    )

    # Plot differentially expressed genes for each cluster
    logger.info("Plot differentially expressed genes for each cluster in elbow plot...")
    sc.pl.rank_genes_groups(
        adata,
        n_genes=10,
        sharey=False,
        ncols=3,
        legend_fontsize=10,
        show=False,
        save="rank_genes_groups_leiden.png",
    )

    # Make a dataframe of marker expression
    logger.info("Save files for differentially expressed genes for each cluster...")
    logger.info("File 1...")
    markers = sc.get.rank_genes_groups_df(adata, None)
    markers = markers[(markers["pvals_adj"] < 0.05) & (markers["logfoldchanges"] > 0.5)]
    markers.to_excel(
        module_dir / "markers.xlsx",
        index=False,
    )

    logger.info("File 2...")
    # Define the number of clusters
    clusters_list = len(adata.obs["leiden"].astype(str).unique())

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

    logger.info("File 3...")
    # Create a dictionary to store DataFrames for each cluster
    cluster_dict = {}
    os.makedirs(os.path.join(module_dir, "cluster_diff_genes"), exist_ok=True)
    for cluster_number in range(clusters_list):
        # print(cluster_number)
        current_cluster = markers[markers["group"] == str(cluster_number)].sort_values(
            by="logfoldchanges", ascending=False
        )  # make a dataframe of the current cluster
        cluster_dict[f"cluster_{cluster_number}"] = (
            current_cluster  # Store the DataFrame in the dictionary
        )
        # Export the DataFrame to a CSV file
        csv_filename = os.path.join(
            output_path,
            module_name,
            "cluster_diff_genes",
            f"cluster_{cluster_number}_data.csv",
        )
        current_cluster.to_csv(csv_filename, index=False)
        print(f"Exported cluster {cluster_number} data to {csv_filename}")

    # Rename the clusters based on the markers
    logger.info("Renaming clusters based on markers...")
    # Get unique clusters
    unique_clusters = (
        adata.obs["leiden"].astype(str).unique()
    )  # Get unique cluster names
    cluster_names = {
        cluster: f"Cluster_{cluster}" for cluster in unique_clusters
    }  # Create a mapping of cluster names
    adata.obs["cell_type"] = (
        adata.obs["leiden"].astype(str).map(cluster_names)
    )  # Map the cluster names to the cell_type column

    # Save anndata object
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")
    logger.info("Annotation module completed successfully.")


if __name__ == "__main__":
    # Set up logger
    configure_logging()
    logger = getLogger("recode_st.3_annotate")  # re-name the logger to match the module

    # Set seed
    seed_everything(21122023)

    try:
        run_annotate()
    except FileNotFoundError as err:
        logger.error(f"File not found: {err}")
