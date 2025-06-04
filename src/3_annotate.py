# Import packages
import warnings  # ? what is the best way to suppress warnings from package inputs?

warnings.filterwarnings("ignore")

import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# import torch
import openpyxl
import pandas as pd
import scanpy as sc
import seaborn as sns
import spatialdata as sd
import squidpy as sq

# from helper_function.py import seed_everything

# Set seed
# seed_everything(21122023)


# Set variables
module_name = "3_annotate"

# Set directories
base_dir = Path(os.getenv("RECODE_BASE_DIR"))
input_path = base_dir
output_path = base_dir / "analysis"
logging_path = output_path / "logging"

# Confirm directories exist
if not Path(input_path).exists():
    raise FileNotFoundError(f"Input path {input_path} does not exist.")
if not Path(output_path).exists():
    raise FileNotFoundError(f"Output path {output_path} does not exist.")

# Create output directories if they do not exist
os.makedirs(output_path / module_name, exist_ok=True)

# Set up logging
os.makedirs(logging_path, exist_ok=True)
logging.basicConfig(
    filename=logging_path / f"{module_name}.txt",  # output file
    filemode="w",  # overwrites the file each time
    format="%(asctime)s - %(levelname)s - %(message)s",  # log format
    level=logging.INFO,  # minimum level to log
)

# Import data
logging.info("Loading Xenium data...")
adata = sc.read_h5ad(Path(output_path) / "2_DR/adata.h5ad")

# change directory to output_path/module_name
os.chdir(
    output_path / module_name
)  # need to so plots save in the correct directory

# Annotate cell clusters
# Calculate the differentially expressed genes for every cluster, compared to the rest of the cells in our adata

logging.info("Calculating differentially expressed genes for each cluster...")
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

logging.info(
    "Plotting the top differentially expressed genes for each cluster..."
)
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby="leiden",
    standard_scale="var",
    n_genes=5,
    show=False,
    save="dotplot.png",
)

# Plot differentially expressed genes for each cluster
logging.info(
    "Plot differentially expressed genes for each cluster in elbow plot..."
)
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
logging.info(
    "Save files for differentially expressed genes for each cluster..."
)
logging.info("File 1...")
markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[
    (markers["pvals_adj"] < 0.05) & (markers["logfoldchanges"] > 0.5)
]
markers.to_excel(
    Path(output_path) / module_name / "markers.xlsx",
    index=False,
)

logging.info("File 2...")
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
    Path(output_path) / module_name / "top_differentially_expressed_genes.csv",
    index=True,
)

logging.info("File 3...")
# Create a dictionary to store DataFrames for each cluster
cluster_dict = {}
os.makedirs(
    os.path.join(output_path, module_name, "cluster_diff_genes"), exist_ok=True
)
for cluster_number in range(clusters_list):
    # print(cluster_number)
    current_cluster = markers[
        markers["group"] == str(cluster_number)
    ].sort_values(
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
logging.info("Renaming clusters based on markers...")
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
adata.write_h5ad(output_path / f"{module_name}/adata.h5ad")
logging.info(f"Data saved to {module_name}/adata.h5ad")
