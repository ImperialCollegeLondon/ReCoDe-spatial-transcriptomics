"""Muspan module."""

import warnings
from logging import getLogger

import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns

from recode_st.config import IOConfig, MuspanModuleConfig
from recode_st.helper_function import configure_scanpy_figures

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)

# TODO: replace "cell_type" with variable


def run_muspan(config: MuspanModuleConfig, io_config: IOConfig):
    """Run Muspan on Xenium data."""
    try:
        import muspan as ms
    except ModuleNotFoundError as err:
        logger.error(
            "Could not load necessary MuSpAn package. You can obtain this with:\n"
            "    pip install 'recode_st[muspan] @ git+"
            "https://github.com/ImperialCollegeLondon/ReCoDe-spatial-transcriptomics.git"
        )
        raise err
    # Set variables
    module_dir = io_config.output_dir / config.module_name
    domain_name = config.domain_name
    transcripts_of_interest = config.transcripts_of_interest
    adata_cell_id = config.adata_cell_id
    cluster_labels = config.cluster_labels

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set figure settings to ensure consistency across all modules
    configure_scanpy_figures(str(io_config.output_dir))
    cmap = sns.color_palette("Spectral", as_cmap=True)

    # Import data
    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(io_config.output_dir / "5_spatial_stats" / "adata.h5ad")

    # Create muspan object
    logger.info("Creating MuSpAn domain object...")
    domain = ms.io.xenium_to_domain(
        path_to_xenium_data=str(io_config.xenium_dir),  # Convert Path to string
        cells_from_selection_csv=str(io_config.area_path),  # Convert Path to string
        domain_name=domain_name,
        load_transcripts=True,
        selected_transcripts=transcripts_of_interest,
        load_nuclei=True,
        load_cells_as_shapes=True,
        exclude_no_nuclei_cells=True,
    )
    logger.info(f"Domain object created: {domain}")

    # Filer adata only to include cells within the area of interest
    logger.info("Filtering adata to include only cells in the area of interest...")
    # Get cell IDs in their original order (preserving duplicates if any)
    domain = map_cell_types_to_domain(adata, domain, adata_cell_id, cluster_labels)

    logger.info("Queries to isolate the different objects within the MuSpAn domain")

    qCells = ms.query.query(
        domain, ("Collection",), "is", "Cell boundaries"
    )  # Query to isolate cell boundaries

    qTrans = ms.query.query(
        domain, ("Collection",), "is", "Transcripts"
    )  # Query to isolate transcripts
    qNuc = ms.query.query(
        domain, ("Collection",), "is", "Nucleus boundaries"
    )  # Query to isolate nucleus boundaries

    logger.info(f"Visualize the MuSpAn domain: {domain_name}")
    fig, ax = plt.subplots(figsize=(20, 15), nrows=2, ncols=2)
    # Visualise all objects in the MuSpAn domain
    ms.visualise.visualise(domain, ax=ax[0, 0], marker_size=0.05)
    ax[0, 0].set_title("All objects")

    ms.visualise.visualise(
        domain,
        color_by=("label", "cell_type"),
        ax=ax[0, 1],
        objects_to_plot=qCells,
    )
    ax[0, 1].set_title("Cell type")

    # Visualise transcripts, colored by 'Transcript'
    ms.visualise.visualise(
        domain,
        color_by=("label", "Transcript ID"),
        ax=ax[1, 0],
        objects_to_plot=qTrans,
        marker_size=5,
    )
    ax[1, 0].set_title("Transcripts")

    # Visualise nuclei, colored by 'Nucleus Area'
    ms.visualise.visualise(
        domain,
        color_by=("label", "Nucleus Area"),
        ax=ax[1, 1],
        objects_to_plot=qNuc,
    )
    ax[1, 1].set_title("Nuclei")

    plt.tight_layout()
    plt.savefig(module_dir / "muspan_domain_visualization.png")
    logger.info("Muspan Domain Visualization plotted and saved")

    logger.info("Convert cell boundaries to cell centres (centroids)")
    domain.convert_objects(
        population=("Collection", "Cell boundaries"),
        object_type="point",
        conversion_method="centroids",
        collection_name="Cell centroids",
        inherit_collections=False,
    )

    # Plot the cell centres with color based on 'cell_type'
    plt.figure(figsize=(10, 6))
    ms.visualise.visualise(
        domain,
        objects_to_plot=("collection", "Cell centroids"),
        color_by="cell_type",
        ax=plt.gca(),
    )
    plt.tight_layout()
    plt.savefig(module_dir / "muspan_cell_centroids.png")
    logger.info(
        "Cell centroids visualized with cell_type color coding plotted and saved"
    )

    # Overlay the cell centres with cell type
    plt.figure(figsize=(10, 6))
    ms.visualise.visualise(
        domain, objects_to_plot=qCells, color_by="cell_type", ax=plt.gca()
    )
    ms.visualise.visualise(
        domain,
        objects_to_plot=("collection", "Cell centroids"),
        color_by=("constant", "black"),
        ax=plt.gca(),
        marker_size=1,
        add_cbar=False,
    )
    plt.tight_layout()
    plt.savefig(module_dir / "muspan_cell_centroids_n_boundaries.png")
    logger.info(
        "Cell centroids & cell boundaries with cell type color plotted and saved"
    )

    # Save domain
    ms.io.save_domain(
        domain,
        name_of_file="muspan_object",
        path_to_save=str(module_dir),
    )
    logger.info("Domain saved")


def map_cell_types_to_domain(adata, domain, adata_cell_id, cluster_labels):
    """Maps cell type (cluster) labels from an AnnData to a domain object on cell ID.

    This function filters the input AnnData object (`adata`) to include only cells
    present in the specified domain, then maps the cluster labels (cell types) from
    `adata` to the domain object, preserving the order and duplicates of cell IDs as
    they appear in the domain. The mapped cell type labels are added to the domain
    object under the specified label name.

    Args:
        adata (AnnData): The annotated data matrix containing cell metadata and
            cluster labels.
        domain: An object representing a spatial or logical domain, with cell IDs
            accessible via `domain.labels["Cell ID"]["labels"]` and a method
            `add_labels` for adding new labels.
        adata_cell_id (str): The column name in `adata.obs` that contains the cell
            IDs.
        cluster_labels (str): The column name in `adata.obs` that contains the
            cluster or cell type labels.

    Returns:
        None. The function modifies the `domain` object in place by adding a
        new label with the mapped cell types.

    Logs:
        - Number of unique and total cell IDs in the domain.
        - Number of cells before and after filtering `adata`.
        - Addition of cell type IDs to the domain.
        - Keys of labels in the domain after addition.
        - Lengths of relevant lists and the count of 'Unknown' cell types.
    """
    domain_cell_ids_ordered = [
        str(cell_id) for cell_id in domain.labels["Cell ID"]["labels"]
    ]

    # Get unique cell IDs for filtering adata
    domain_cell_ids_unique = set(domain_cell_ids_ordered)

    logger.info(f"Number of unique cells in the domain: {len(domain_cell_ids_unique)}")
    logger.info(
        f"Total cell entries in domain (including duplicates): "
        f"{len(domain_cell_ids_ordered)}"
    )

    # Filter adata to include only cells in the area of interest
    filt_adata = adata[
        adata.obs[adata_cell_id].isin(domain_cell_ids_unique)
    ]  # filter adata object

    logger.info(f"Filtered adata from {adata.n_obs} to {filt_adata.n_obs} cells")

    # Add cell cluster IDs
    logger.info("Adding cell_type IDs to domain with cluster labels")

    # Create a mapping from cell_id to cell_type
    cell_id_to_type = dict(
        zip(filt_adata.obs[adata_cell_id], filt_adata.obs[cluster_labels])
    )

    # Get cell types in the same order as domain cell IDs
    # (preserving order and duplicates)
    cell_types_ordered = [
        cell_id_to_type.get(cell_id, "Unknown") for cell_id in domain_cell_ids_ordered
    ]

    # Add cell_type label to the domain
    domain.add_labels(label_name=cluster_labels, labels=cell_types_ordered)
    logger.info(f"Label keys in domain: {domain.labels.keys()}")
    logger.info(f"Length of cell_types_ordered: {len(cell_types_ordered)}")
    logger.info(f"Length of domain cell IDs: {len(domain_cell_ids_ordered)}")
    logger.info(
        f"Number of 'Unknown' cell types: {cell_types_ordered.count('Unknown')}"
    )

    # Return muspan domain with cell types mapped
    return domain
