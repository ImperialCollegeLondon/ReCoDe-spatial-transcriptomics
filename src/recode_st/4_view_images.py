"""Image viewing module."""

# Import packages
import logging
import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

from recode_st.paths import logging_path, output_path

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    # Set variables
    module_name = "4_view_images"  # name of the module
    gene_list = ["EPCAM", "CD3D", "CD68", "VWF", "PTPRC", "ACTA2"]
    module_dir = output_path / module_name

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set up logging
    os.makedirs(
        logging_path, exist_ok=True
    )  # should set up all these directories at the start of the pipeline?
    logging.basicConfig(
        filename=Path(logging_path) / f"{module_name}.txt",  # output file
        filemode="w",  # overwrites the file each time
        format="%(asctime)s - %(levelname)s - %(message)s",  # log format
        level=logging.INFO,  # minimum level to log
    )

    # change directory to output_path/module_name
    os.chdir(module_dir)
    logging.info(f"Changed directory to {module_dir}")

    # Import data
    logging.info("Loading Xenium data...")
    adata = sc.read_h5ad(module_dir / "3_annotate/adata.h5ad")

    # View plots
    logging.info("Visualize clusters on tissue...")
    sq.pl.spatial_scatter(
        adata,
        library_id="spatial",
        shape=None,
        outline=False,
        color=["leiden", "total_counts"],
        wspace=0.4,
        size=1,
    )
    plt.tight_layout()
    plt.savefig(module_dir / f"{module_name}/images.png", dpi=300)
    plt.close()
    logging.info(f"Saved plots to {module_dir / 'images.png'}")

    # View specific gene expression
    logging.info("Plotting genes of interest on tissue...")
    sq.pl.spatial_scatter(
        adata,
        library_id="spatial",
        color=gene_list,
        shape=None,
        size=2,
        img=False,
        save="gene_expression.png",
    )

    # Save anndata object
    adata.write_h5ad(module_dir / f"{module_name}/adata.h5ad")
    logging.info(f"Data saved to {module_dir / 'adata.h5ad'}")
    logging.info("Imaging module completed successfully.")
