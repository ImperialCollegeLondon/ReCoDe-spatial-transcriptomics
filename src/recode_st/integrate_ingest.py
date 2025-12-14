"""Integrate scRNAseq and STx using ingest."""

import warnings
from logging import getLogger
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from recode_st.config import IntegrateIngestModuleConfig, IOConfig
from recode_st.helper_function import configure_scanpy_figures

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)

# Define global variables - SHOULD THESE BE HERE OR INSIDE THE RUN FUNCTION?
INGEST_LABEL_COL = "ingest_pred_cell_type"
REF_CELL_LABEL_COL = "cell_type"  # Column in reference data with cell type labels
BATCH_COL = "dataset_origin"
ANNOTATION_COLS = [  # Define the annotation columns to transfer
    "transf_ann_level_1_label",
    "transf_ann_level_2_label",
    "transf_ann_level_3_label",
    "transf_ann_level_4_label",
    "transf_ann_level_5_label",
]
HLCA_INT_SAVE = Path(
    "/rds/general/user/sep22/ephemeral/recode_hlca_full_processed.h5ad"
)
SUBSET_HLCA_INT_SAVE = Path(
    "/rds/general/user/sep22/ephemeral/recode_subset_hlca_full_processed.h5ad"
)


def prepare_integrated_datasets(gene_id_dict_path, adata_ref, adata):
    """Ensures ST and reference dataset contain are in correct format for integration.

    Both datasets should have:
    1. The same list of genes
    2. Labeled by the ensembl ID rather than gene symbol
    3. Genes are in the same order.
    This needs to be done prior to integration!

    Args:
        gene_id_dict_path (str | Path): Path to gene ID dictionary .csv file.
        First column: named "gene_symbol" containing the gene ID for all genes in
        STx dataset.
        Second column:  named "ensembl_ID" containing the ensembl ID of the
        corresponding gene.
        adata_ref (anndata.AnnData)): Reference scRNAseq AnnData object.
        adata (anndata.AnnData): STx AnnData object.

    Raises:
        FileNotFoundError: _description_
        ValueError: _description_

    Returns:
        anndata.AnnData: formatted adata_ref_subset and adata_ingest AnnData objects.

    """
    logger.info("Confirm Xenium data and reference data have the same genes...")
    # Replace ensembl ID with gene symbols from adata_ref for matching
    try:
        gene_id_dict = pd.read_csv(
            gene_id_dict_path, index_col=0
        )  # dictionary with ensembl and gene symbols
    except FileNotFoundError as e:
        logger.error(f"Gene ID dictionary file not found: {gene_id_dict_path}")
        raise FileNotFoundError(
            f"Gene ID dictionary file not found: {gene_id_dict_path}"
        ) from e
    except pd.errors.ParserError as e:
        logger.error(f"Error parsing gene ID dictionary file: {gene_id_dict_path}")
        raise ValueError(
            f"Error parsing gene ID dictionary file: {gene_id_dict_path}"
        ) from e
    except Exception as e:
        logger.error(f"Error reading gene ID dictionary file: {gene_id_dict_path}: {e}")
        raise

    # Add ensembl_id to STx data
    adata.var["ensembl_id"] = adata.var.index.map(gene_id_dict["ensembl_id"])
    missing = adata.var["ensembl_id"].isna().sum()
    if missing > 0:
        logger.warning(f"{missing} genes in STx have no matching Ensembl ID.")

    # List of genes shared between datasets based on ensembl IDs
    var_names = adata_ref.var_names.intersection(adata.var["ensembl_id"])
    logger.info(f"Number of common genes: {len(var_names)}")
    if len(var_names) < 500:
        logger.warning(
            "Warning: Less than 500 common genes found between datasets. "
            "Integration may not perform well."
        )

    # Subset STx data to common genes
    mask = adata.var["ensembl_id"].isin(
        var_names
    )  # Create mask to filter genes based on ensembl IDs
    adata_ingest = adata[:, mask].copy()

    adata_ingest.var_names = adata_ingest.var[
        "ensembl_id"
    ]  # Rename var_names to the Ensembl IDs

    # Check for duplicate gene names in ST dataset
    if adata_ingest.var_names.has_duplicates:
        raise ValueError("Duplicate Ensembl IDs detected in ST dataset after mapping.")

    logger.info(f"First 5 gene names in ST: {adata_ingest.var_names[:5].tolist()}")
    # Subset reference datasets to common genes
    adata_ref_subset = adata_ref[:, var_names].copy()

    # After subsetting both datasets to common genes, add these checks:
    # 1. Check the sets of genes are the same
    ref_genes = set(adata_ref_subset.var_names)
    st_genes = set(adata_ingest.var_names)
    same_set = ref_genes == st_genes
    logger.info(f"Gene sets identical: {same_set}")

    if not same_set:
        missing_in_ref = st_genes - ref_genes
        missing_in_ingest = ref_genes - st_genes
        logger.warning(f"Genes in ST but not in reference: {len(missing_in_ref)}")
        logger.warning(f"Genes in reference but not in ST: {len(missing_in_ingest)}")
        raise ValueError("Datasets do not have identical gene sets. Cannot reorder!")

    # 2. Check order
    same_order = adata_ref_subset.var_names.equals(adata_ingest.var_names)
    logger.info(f"Genes in same order: {same_order}")

    # 3. If order doesn't match, reorder adata_ingest to match adata_ref
    if not same_order:
        logger.info("Reordering adata_ingest genes to match reference...")
        adata_ingest = adata_ingest[:, adata_ref_subset.var_names].copy()
        logger.info("Genes reordered successfully.")

    # Confirm shapes
    logger.info(f"Reference dataset shape: {adata_ref_subset.shape}")
    logger.info(f"ST dataset shape: {adata_ingest.shape}")

    return adata_ref_subset, adata_ingest


def process_reference_data(config, io_config, adata_ref):
    """Preprocesses the reference scRNA-seq dataset if PCA and UMAP are missing.

    This function checks whether the reference AnnData object already processed.
    If not, it performs standard single-cell preprocessing steps, including
    normalization, log-transformation, highly variable gene selection, PCA,
    neighbor graph construction, and UMAP embedding.
    The resulting AnnData object is saved to disk and optionally visualized.

    Args:
        config (SimpleNamespace or dict): Configuration object containing
            module-specific parameters, including the module name used to
            construct figure filenames.
        io_config (SimpleNamespace or dict): I/O configuration object with paths
            to input and output data. Must include ``ref_path`` to save the
            processed reference dataset.
        adata_ref (anndata.AnnData): Reference single-cell RNA-seq dataset (e.g.,
            HLCA) to be processed or verified. The object is modified in place.

    Side Effects:
        - Writes the processed reference AnnData object to `HLCA_INT_SAVE`.
        - Saves a UMAP visualization as a PNG file.

    Logs:
        - Whether preprocessing is needed.
        - Each major preprocessing step.
        - Path to the saved output file.

    Returns:
        adata_ref (anndata.AnnData): Processed reference AnnData object.
    """
    logger.info("Checking if reference scRNA-seq data needs processing...")

    # Check if we need to preprocess the reference data
    if HLCA_INT_SAVE.exists():
        logger.info(f"Processed HLCA reference data found at {HLCA_INT_SAVE}.")
        logger.info(f"Loading existing processed data from {HLCA_INT_SAVE}.")
        adata_ref = sc.read_h5ad(HLCA_INT_SAVE)
    else:
        logger.info("Preprocessing HLCA reference data...")

        # Basic filtering and normalization
        sc.pp.filter_cells(adata_ref, min_genes=200)
        sc.pp.filter_genes(adata_ref, min_cells=10)

        # Checking for counts layer
        if "counts" in adata_ref.layers:
            logger.info("Counts layer present!")
        else:
            logger.warning("Counts layer not found in adata_ref.layers!")

        # Check if data is already normalized
        logger.info("Checking normalization...")
        totals = adata_ref.X.sum(axis=1)
        totals = totals.A1 if hasattr(totals, "A1") else np.array(totals).ravel()

        if np.median(totals) > 2e4:  # crude but works for most scRNA-seq
            logger.info("Data does not appear to be normalized! Normalizing...")
            sc.pp.normalize_total(adata_ref, target_sum=1e4)
            logger.info("Normalized total counts.")
        else:
            logger.info("Data already appears normalized. Skipping.")

        # Check if data is log-transformed
        logger.info("Checking log transformation...")
        X_sample = adata_ref.X[:100, :100]
        arr = X_sample.A if hasattr(X_sample, "A") else X_sample

        if np.allclose(arr, np.round(arr)):  # looks like raw counts
            logger.info(
                "Data does not appear to be log transformed! Transforming data..."
            )
            sc.pp.log1p(adata_ref)
            logger.info("Applied log1p transformation.")
        else:
            logger.info("Data already appears log-transformed. Skipping.")

        logger.info("Calculating HVG...")
        # Feature selection and scaling
        sc.pp.highly_variable_genes(adata_ref, n_top_genes=5000, flavor="seurat_v3")
        # sc.pp.scale(adata_ref, max_value=10) - scaling skipped for large datasets

        # Dimensionality reduction
        logger.info("Computing PCA...")
        sc.tl.pca(adata_ref, n_comps=75, svd_solver="arpack")
        sc.pp.neighbors(adata_ref, n_neighbors=30, n_pcs=75)
        logger.info("Computing UMAP...")
        sc.tl.umap(adata_ref)

        # Save processed reference
        adata_ref.write_h5ad(HLCA_INT_SAVE)
        logger.info(
            f"Finished preprocessing. Processed HLCA reference saved to {HLCA_INT_SAVE}"
        )

    logger.info("Visualizing reference data UMAP...")
    # Plot UMAP
    sc.pl.umap(
        adata_ref,
        color=REF_CELL_LABEL_COL,
        title="HLCA reference data UMAP",
        save=f"_{config.module_name}_hlca_umap.png",
    )
    return adata_ref


def process_subset_reference_data(config, io_config, adata_ref_subset):
    """Preprocesses the reference scRNA-seq dataset if PCA and UMAP are missing.

    This function checks whether the reference AnnData object already processed.
    If not, it performs standard single-cell preprocessing steps, including
    normalization, log-transformation, highly variable gene selection, PCA,
    neighbor graph construction, and UMAP embedding.
    The resulting AnnData object is saved to disk and optionally visualized.

    Args:
        config (SimpleNamespace or dict): Configuration object containing
            module-specific parameters, including the module name used to
            construct figure filenames.
        io_config (SimpleNamespace or dict): I/O configuration object with paths
            to input and output data. Must include ``ref_path`` to save the
            processed reference dataset.
        adata_ref_subset (anndata.AnnData): Subset of reference single-cell RNA-seq
        dataset (e.g., HLCA) to be processed or verified.
        The object is modified in place.

    Side Effects:
        - Writes the processed reference AnnData object to `SUBSET_HLCA_INT_SAVE`.
        - Saves a UMAP visualization as a PNG file.

    Logs:
        - Whether preprocessing is needed.
        - Each major preprocessing step.
        - Path to the saved output file.

    Returns:
        adata_ref_subset (anndata.AnnData): Processed reference AnnData object.
    """
    logger.info("Checking if reference scRNA-seq data needs processing...")

    # Check if we need to preprocess the reference data
    if SUBSET_HLCA_INT_SAVE.exists():
        logger.info(f"Processed HLCA reference data found at {SUBSET_HLCA_INT_SAVE}.")
        logger.info(f"Loading existing processed data from {SUBSET_HLCA_INT_SAVE}.")
        adata_ref_subset = sc.read_h5ad(SUBSET_HLCA_INT_SAVE)
    else:
        logger.info("Preprocessing subset of HLCA reference data...")
        logger.info("Compute PCA and neighbors for ingest reference subset...")
        sc.tl.pca(adata_ref_subset, n_comps=75, svd_solver="arpack")
        sc.pp.neighbors(adata_ref_subset, n_neighbors=30, n_pcs=75)
        sc.tl.umap(adata_ref_subset)
        sc.pl.umap(
            adata_ref_subset,
            color=REF_CELL_LABEL_COL,
            title="HLCA reference subset for ingest UMAP",
            save=f"_{config.module_name}_hlca_subset_ingest_umap.png",
        )
        logger.info(f"UMAP for {REF_CELL_LABEL_COL} saved for adata_ref_subset")

        # Save processed reference
        adata_ref_subset.write_h5ad(SUBSET_HLCA_INT_SAVE)
        logger.info("Finished preprocessing")
        logger.info(f"Processed HLCA reference saved to {SUBSET_HLCA_INT_SAVE}")
    return adata_ref_subset


def ingest_integration(adata_ref, adata, adata_ingest):
    """Integrates STx data with a reference scRNA-seq dataset using ingest.

    This function performs label transfer from a reference single-cell dataset
    (e.g., HLCA) to a STx dataset (e.g., Xenium).
    The ingest approach uses Scanpys `tl.ingest` to
    project query data into the references PCA space and transfer cell-type labels.

    Args:
        adata_ref (anndata.AnnData): Reference scRNA-seq AnnData object containing
            precomputed embeddings (PCA or UMAP) and cell-type annotations.
        adata (anndata.AnnData): Original STx dataset to which
            transferred labels will be written.
        adata_ingest (anndata.AnnData): Copy of the spatial dataset formatted for
            ingestion, aligned by gene set with the reference.

    Returns:
        adata (anndata.AnnData): Original STx dataset with transferred labels added.
    """
    logger.info("Starting integration...")
    logger.info("Integrating data using ingest...")

    # Run ingest to map Xenium data onto HLCA reference
    sc.tl.ingest(
        adata_ingest,
        adata_ref,
        obs=REF_CELL_LABEL_COL,  # Annotation column to use in adata_ref.obs
        # For core HLCA, "cell_type"
        embedding_method="pca",  # or 'umap'
        labeling_method="knn",  # how to transfer labels
        neighbors_key=None,  # use default neighbors from reference
        inplace=True,
    )

    # Rename predicted cell type column
    adata_ingest.obs[INGEST_LABEL_COL] = adata_ingest.obs[REF_CELL_LABEL_COL]
    del adata_ingest.obs[REF_CELL_LABEL_COL]

    # Copy cell type predictions back to original adata
    adata.obs[INGEST_LABEL_COL] = adata_ingest.obs.loc[
        adata.obs_names, INGEST_LABEL_COL
    ]
    logger.info("Ingest integration complete.")

    return adata


def transfer_hierarchical_annotations(
    adata_ref, adata, ANNOTATION_COLS=ANNOTATION_COLS
):
    """Transfer hierarchical annotation levels from reference to spatial data.

    This function creates a mapping from cell types to their hierarchical annotations
    and transfers them to the spatial dataset based on the predicted cell types.

    Args:
        adata_ref (anndata.AnnData): Reference dataset containing the hierarchical
            annotation columns.
        adata (anndata.AnnData): Spatial dataset with predicted cell types that will
            receive the transferred annotations.
        ANNOTATION_COLS (list of str): List of hierarchical annotation column names
            in the reference dataset to be transferred.

    Returns:
        adata (anndata.AnnData): Spatial dataset with added hierarchical annotations.
    """
    logger.info(
        "Transferring hierarchical annotations from reference to spatial data..."
    )

    # Check which columns exist in reference data
    available_cols = [col for col in ANNOTATION_COLS if col in adata_ref.obs.columns]
    missing_cols = [col for col in ANNOTATION_COLS if col not in adata_ref.obs.columns]

    if missing_cols:
        logger.warning(f"Missing annotation columns in reference: {missing_cols}")

    if not available_cols:
        logger.error("No hierarchical annotation columns found in reference data!")
        return adata

    logger.info(f"Found annotation columns: {available_cols}")

    # Check if spatial data has the predicted cell type column
    if INGEST_LABEL_COL not in adata.obs.columns:
        logger.error(
            f"Spatial data missing {INGEST_LABEL_COL} column! Run integration first."
        )
        return adata

    # Create mapping from cell type to hierarchical annotations
    # Group by cell type and take the most common annotation for each level
    logger.info("Creating cell type to annotation mappings...")

    mappings = {}
    for col in available_cols:
        # For each cell type, find the most common annotation at each level
        cell_type_mapping = {}

        for cell_type in adata_ref.obs[REF_CELL_LABEL_COL].unique():
            mask = adata_ref.obs[REF_CELL_LABEL_COL] == cell_type
            annotations = adata_ref.obs.loc[mask, col]

            # Get the most common annotation for this cell type
            if len(annotations) > 0:
                most_common = annotations.mode()
                if len(most_common) > 0:
                    cell_type_mapping[cell_type] = most_common.iloc[0]
                else:
                    cell_type_mapping[cell_type] = "Unknown"
            else:
                cell_type_mapping[cell_type] = "Unknown"

        mappings[col] = cell_type_mapping
        logger.info(f"Created mapping for {col}: {len(cell_type_mapping)} cell types")

    # Apply mappings to spatial data
    logger.info("Applying mappings to spatial data...")

    for col in available_cols:
        # Map predicted cell types to hierarchical annotations
        adata.obs[col] = adata.obs[INGEST_LABEL_COL].map(mappings[col])

        # Fill any unmapped values with "Unknown"
        adata.obs[col] = adata.obs[col].fillna("Unknown")

        # Convert to categorical for efficiency
        adata.obs[col] = adata.obs[col].astype("category")

        # Log statistics
        n_mapped = (adata.obs[col] != "Unknown").sum()
        n_total = len(adata.obs[col])
        logger.info(
            f"Transferred {col}: {n_mapped}/{n_total} cells mapped successfully"
        )

    # Log final results
    logger.info("Hierarchical annotation transfer completed!")
    logger.info(f"Added columns to spatial data: {available_cols}")

    # Show example of the mapping
    if len(available_cols) > 0:
        sample_data = adata.obs[[INGEST_LABEL_COL, *available_cols]].head(3)
        logger.info(f"Sample annotations:\n{sample_data.to_string()}")

    return adata


def visualize_integration(config, cmap, adata):
    """Visualize integration results using UMAPs.

    Args:
        config (SimpleNamespace or dict): Configuration object containing
            module-specific parameters, including the module name used to
            construct figure filenames.
        cmap (_type_): Color map for visualization.
        adata (anndata.AnnData): Original STx dataset to which
            transferred labels will be written.
    """
    logger.info("Generating UMAPs...")
    logger.info(f"Columns in adata: {list(adata.obs.columns)}")

    # Check if UMAP exists, if not compute it
    if "X_umap" not in adata.obsm:
        logger.info("Computing UMAP for visualization...")
        if "X_pca" not in adata.obsm:
            logger.info("Computing PCA first...")
            sc.tl.pca(adata, n_comps=50)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
        sc.tl.umap(adata)

    for method in [INGEST_LABEL_COL]:
        if method not in adata.obs.columns:
            print(f"Missing column: {method}")
            continue

        sc.pl.umap(
            adata,
            color=method,
            title=f"HLCA integration: {method}",
            show=False,  # change to True if you want inline display
            cmap=cmap,
            save=f"_{config.module_name}_{method}.png",
        )

    for obs in ANNOTATION_COLS:
        if obs not in adata.obs.columns:
            print(f"Missing column: {obs}")
            continue

        sc.pl.umap(
            adata,
            color=obs,
            title=f"HLCA integration: {obs}",
            show=False,  # change to True if you want inline display
            cmap=cmap,
            save=f"_{config.module_name}_{obs}.png",
        )

    logger.info(f"UMAP plots saved to {sc.settings.figdir}")


def run_integration(config: IntegrateIngestModuleConfig, io_config: IOConfig):
    """Integrate scRNAseq and STx data using scANVI and ingest.

    adata_ref_subset and adata_ingest are used for ingest integration.
    adata_ref and data are used for scANVI integration.

    Args:
        config (IntegrateModuleConfig): Integration module configuration object.
        io_config (IOConfig): IO configuration object.

    Returns:
        None
    """
    # Variables

    # Name of the column to store label transfer results in adata.obs
    module_dir = io_config.output_dir / config.module_name

    # Paths to input data
    ref_path = io_config.ref_path
    gene_id_dict_path = io_config.gene_id_dict_path

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set figure directory for this module (overrides global setting)
    sc.settings.figdir = module_dir

    # Set figure settings to ensure consistency across all modules
    configure_scanpy_figures(str(io_config.output_dir))
    cmap = sns.color_palette("Spectral", as_cmap=True)

    logger.info("Starting integration of scRNAseq and spatial transcriptomics data...")

    logger.info("Loading scRNAseq data from HLCA ...")
    adata_ref = sc.read_h5ad(ref_path)

    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(io_config.output_dir / "2_dimension_reduction" / "adata.h5ad")

    logger.info("Processing reference data...")
    process_reference_data(config, io_config, adata_ref)

    # 1. INTEGRATION USING INGEST
    logger.info("Formatting data for ingest integration...")
    adata_ref_subset, adata_ingest = prepare_integrated_datasets(
        gene_id_dict_path, adata_ref, adata
    )

    logger.info("Processing subset of reference data for ingest...")
    adata_ref_subset = process_subset_reference_data(
        config, io_config, adata_ref_subset
    )

    logger.info("Starting integration methods...")

    # 1. Perform ingest integration
    logger.info("Performing integration using: sc.tl.ingest...")
    ingest_integration(adata_ref_subset, adata, adata_ingest)

    # 2. Transfer hierarchical annotations
    logger.info("Transferring hierarchical annotations...")
    transfer_hierarchical_annotations(adata_ref, adata)

    # Delete adata_ref_subset and adata_ingest to free memory
    logger.info("Delete adata_ref_subset and adata_ingest to free memory...")
    del adata_ref_subset
    del adata_ingest

    logger.info("Visualizing integration results...")
    visualize_integration(config, cmap, adata)

    logger.info("Saving integrated data...")
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Integrated data saved to {module_dir}.")

    logger.info("Ingest integration complete.")
