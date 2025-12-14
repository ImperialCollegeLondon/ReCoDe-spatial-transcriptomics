"""Integrate scRNAseq and STx."""

import gc
import warnings
from logging import getLogger
from pathlib import Path

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scvi.model import SCANVI, SCVI

from recode_st.config import Config, IntegrateSCVIModuleConfig, IOConfig

# Suppress specific warnings to reduce noise in logs
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)

# Define global constants for integration
INGEST_LABEL_COL = "ingest_pred_cell_type"
SCANVI_LABEL_COL = "scANVI_pred_cell_type"
REF_CELL_LABEL_COL = "cell_type"  # Column in reference data with cell type labels
BATCH_COL = "dataset_origin"
REFERENCE_KEY = "Ref"
SPATIAL_KEY = "STx"
SCVI_LATENT_KEY = "X_scVI"
SCANVI_LATENT_KEY = "X_scANVI"
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


def verify_counts_layer(adata, data_type="data"):
    """Verify that counts layer exists and contains valid count data.

    Args:
        adata (anndata.AnnData): AnnData object to verify
        data_type (str): Description of data type for logging
            (e.g., "reference", "spatial")

    Raises:
        ValueError: If counts layer is missing or invalid
    """
    logger.info(f"Verifying {data_type} data integrity...")

    if "counts" in adata.layers:
        logger.info(f"✓ {data_type.capitalize()} data has counts layer")
        logger.info(f"✓ {data_type.capitalize()} layers: {list(adata.layers.keys())}")
        logger.info(f"✓ Counts layer shape: {adata.layers['counts'].shape}")
        logger.info(f"✓ {data_type.capitalize()} data shape: {adata.shape}")

        # Verify counts layer matches main data dimensions
        if adata.layers["counts"].shape != adata.X.shape:
            counts_shape = adata.layers["counts"].shape
            x_shape = adata.X.shape
            logger.error(
                f"✗ Counts layer shape mismatch: counts={counts_shape}, X={x_shape}"
            )
            raise ValueError("Counts layer dimensions don't match main data matrix")

        # Validate counts data quality
        logger.info("Validating counts data quality...")
        counts_sample = adata.layers["counts"][:100, :100]
        counts_arr = counts_sample.A if hasattr(counts_sample, "A") else counts_sample
        is_integer = np.allclose(counts_arr, np.round(counts_arr))
        median_total = np.median(counts_arr.sum(axis=1))
        has_positive_values = (counts_arr >= 0).all()

        if is_integer and median_total > 100 and has_positive_values:
            logger.info("✓ Counts layer contains valid integer count data")
            logger.info(f"✓ Median counts per cell: {median_total:.0f}")
        else:
            logger.warning(
                f"⚠ Counts layer validation issues: integer-like={is_integer}, "
                f"median_total={median_total:.1f}, non_negative={has_positive_values}"
            )
            # Don't fail here, but warn user
            if not has_positive_values:
                logger.error("✗ Counts layer contains negative values!")
                raise ValueError("Invalid counts layer: contains negative values")
    else:
        logger.error(f"✗ {data_type.capitalize()} data missing 'counts' layer!")
        logger.error("Available layers: " + str(list(adata.layers.keys())))
        if data_type.lower() == "reference":
            logger.error(
                "Please run hlca_full_filt_process.py first to create "
                "processed data with counts layer"
            )
        else:
            logger.error(
                "Please ensure your data preprocessing preserves the counts layer"
            )
        raise ValueError(
            f"{data_type.capitalize()} data missing required 'counts' layer"
        )


def subset_reference(
    base_config,
    adata_ref: anndata.AnnData,
    target_cells: int = 200000,
    stratify_col: str = REF_CELL_LABEL_COL,
    min_per_type: int = 100,
    sampling_strategy="proportional",
    rare_cell_boost: float = 2.0,
):
    """Subsample reference data while preserving cell type diversity.

    Args:
        base_config (SimpleNamespace or dict): Configuration object used to set seed.
        adata_ref (anndata.AnnData): Reference AnnData object.
        target_cells (int): Target number of cells after subsampling.
        stratify_col (str): Column in adata_ref.obs to stratify sampling by.
        min_per_type (int): Minimum cells to sample per cell type.
        sampling_strategy (str): One of:
            - "proportional": Sample proportionally to original frequencies
            - "balanced": Give equal weight to all cell types
            - "sqrt": Compromise between proportional and balanced (sqrt of proportions)
        rare_cell_boost (float): Multiplier for min_per_type for rare cell types.

    Returns:
        anndata.AnnData: Subsampled reference AnnData object.
    """
    rng = np.random.default_rng(base_config.seed)  # set just in case
    cell_type_counts = adata_ref.obs[stratify_col].value_counts()
    n_types = len(cell_type_counts)  # number of unique cell types
    total_cells = len(adata_ref)  # total number of cells in reference

    # Calculate sampling targets based on strategy
    if sampling_strategy == "proportional":
        # Keep original proportions
        sampling_fracs = cell_type_counts / total_cells
        base_targets = (sampling_fracs * target_cells).round().astype(int)

    elif sampling_strategy == "balanced":
        # Equal representation for all cell types
        base_targets = pd.Series(target_cells // n_types, index=cell_type_counts.index)

    elif sampling_strategy == "sqrt":
        # Square root transformation - compromise between proportional and balanced
        sqrt_counts = np.sqrt(cell_type_counts)
        sampling_weights = sqrt_counts / sqrt_counts.sum()
        base_targets = (sampling_weights * target_cells).round().astype(int)

    else:
        raise ValueError(f"Unknown sampling_strategy: {sampling_strategy}")

    # Apply minimum per type, with boost for very rare cell types
    rare_threshold = total_cells * 0.001  # <0.1% of total
    final_targets = {}

    for ct, count in cell_type_counts.items():
        base_target = base_targets[ct]

        # Determine minimum based on rarity
        if count < rare_threshold:
            effective_min = int(min_per_type * rare_cell_boost)
        else:
            effective_min = min_per_type

        # Don't sample more than available, but ensure minimum
        final_targets[ct] = min(count, max(base_target, effective_min))

    # Adjust if we're over target_cells due to minimums
    total_target = sum(final_targets.values())
    if total_target > target_cells:
        # Scale down proportionally, but preserve minimums
        excess = total_target - target_cells
        adjustable_cts = {
            ct: target for ct, target in final_targets.items() if target > min_per_type
        }

        if adjustable_cts:
            # Remove excess proportionally from adjustable cell types
            total_adjustable = sum(adjustable_cts.values())
            for ct in adjustable_cts:
                reduction = int(excess * (adjustable_cts[ct] / total_adjustable))
                final_targets[ct] = max(min_per_type, final_targets[ct] - reduction)

    # Perform sampling
    sampled_indices = []
    for ct, n_sample in final_targets.items():
        idx = adata_ref.obs.loc[adata_ref.obs[stratify_col] == ct].index
        sampled = rng.choice(idx, size=n_sample, replace=False)
        sampled_indices.extend(sampled)

    # Shuffle to avoid any ordering bias
    rng.shuffle(sampled_indices)

    adata_subsampled = adata_ref[sampled_indices].copy()

    # Report sampling statistics
    print(f"Subsampled from {total_cells:,} to {len(adata_subsampled):,} cells")
    print(f"Cell types: {n_types}")
    print(f"Sampling strategy: {sampling_strategy}")
    print("\nCell type distribution:")

    orig_pcts = (cell_type_counts / total_cells * 100).round(2)
    new_counts = adata_subsampled.obs[stratify_col].value_counts()
    new_pcts = (new_counts / len(adata_subsampled) * 100).round(2)

    for ct in cell_type_counts.index[:10]:  # Show top 10
        print(
            f"  {ct}: {cell_type_counts[ct]:,} ({orig_pcts[ct]:.2f}%) → "
            f"{new_counts.get(ct, 0):,} ({new_pcts.get(ct, 0):.2f}%)"
        )

    if n_types > 10:
        print(f"  ... and {n_types - 10} more cell types")

    return adata_subsampled


def scVI_integration_check(adata, batch_key=BATCH_COL, cell_type=REF_CELL_LABEL_COL):
    """Check dataset readiness for scVI/scANVI integration.

    Args:
        adata (anndata.AnnData): AnnData object to validate.
        batch_key (str): Column name in adata.obs for batch ID.
        cell_type (str): Column name in adata.obs for cell type labels.

    Returns:
        bool: True if dataset is ready, False otherwise.
    """
    issues = []

    # Expression matrix
    if adata.X is None:
        issues.append("No expression matrix (adata.X) present")

    # Raw count layer
    if "counts" not in adata.layers:
        issues.append("Missing 'counts' layer in adata.layers (required for scVI)")

    # Batch
    if batch_key not in adata.obs.columns:
        issues.append(f"Missing batch column '{batch_key}' in adata.obs")

    if cell_type not in adata.obs.columns:
        issues.append(f"Missing cell type column '{cell_type}' in adata.obs")

    # Outcome
    if issues:
        print("NOT READY FOR scANVI:")
        for issue in issues:
            print(f"  - {issue}")
        return False

    logger.info("Dataset is ready for scVI/scANVI integration.")
    return True


def scVI_integration(config, adata_combined, module_dir):
    """Integrates STx data with a reference scRNA-seq dataset using scVI.

    This function harmonizes a reference single-cell dataset
    (e.g., HLCA) with a STx dataset (e.g., Xenium).

    Args:
        config (SimpleNamespace or dict): Configuration object containing
            module-specific parameters, including the module name used to
            construct figure filenames.
        adata_combined (anndata.AnnData): Combined AnnData object of reference
        and query (STx) data.
        module_dir (Path): Output directory for saving results

    Returns:
        tuple: (adata_combined, scanvi_model) - Combined dataset and trained model
    """
    logger.info("Integrating data using scVI...")

    # Training parameters
    MAX_EPOCHS_SCVI = 200

    # Setup scVI
    logger.info("Setting up scVI model...")
    SCVI.setup_anndata(
        adata_combined,
        layer="counts",  # must be counts layer
        batch_key=BATCH_COL,  # variable you want to perform harmonization over
    )

    logger.info("Initializing scVI model...")
    scvi_model = SCVI(adata_combined, n_layers=3, n_latent=50, n_hidden=128)
    logger.info("Training SCVI model...")
    scvi_model.train(
        max_epochs=MAX_EPOCHS_SCVI,
        batch_size=1024,  # Increase from 128
        plan_kwargs={"lr": 1e-3},  # Explicit learning rate
        early_stopping=True,
        early_stopping_patience=10,
    )

    logger.info("Obtain and visualize latent representation...")
    adata_combined.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation()
    sc.pp.pca(adata_combined, n_comps=50)
    sc.pp.neighbors(adata_combined, use_rep=SCVI_LATENT_KEY, n_neighbors=15)
    sc.tl.umap(adata_combined, min_dist=0.3, spread=1.0)
    sc.pl.umap(
        adata_combined,
        color=[BATCH_COL, REF_CELL_LABEL_COL],
        frameon=False,
        ncols=1,
        save=f"_{config.module_name}_scvi_umap.png",
    )

    logger.info("Saving scANVI model...")
    scvi_model.save(module_dir / "_scvi_ref", overwrite=True)

    logger.info("scVI model saved and available for downstream analysis")
    epochs_completed = len(scvi_model.history["elbo_train"])
    logger.info(f"Model training history: {epochs_completed} epochs completed")

    return adata_combined, scvi_model


def scANVI_label_transfer(config, adata_combined, scvi_model, module_dir):
    """Integrates STx data with a reference scRNA-seq dataset using scVI.

    This function performs label transfer from a reference single-cell dataset
    (e.g., HLCA) to a STx dataset (e.g., Xenium).
    scANVI uses a semi-supervised approach to transfer cell-type labels.

    Args:
        adata_combined (anndata.AnnData): Reference scRNA-seq AnnData object containing
            precomputed embeddings (PCA or UMAP) and cell-type annotations.
            Should have matching genes with adata_ingest.
        scvi_model (scvi.model.SCVI): Trained scVI model on combined data.
        config (SimpleNamespace or dict): Configuration object containing
        module_dir (Path): Output directory for saving results

    Returns:
        tuple: (adata_combined, scanvi_model) - Combined dataset and trained model
    """
    # Set constant variables
    MAX_EPOCHS_SCANVI = 200

    # Format labels for scANVI
    # Create scANVI labels (reference has labels, spatial is 'Unknown')
    SCANVI_CELLTYPE_KEY = (
        "training_celltype_scanvi"  # new column that scANVI will use training labels.
    )
    UNLABELED_CATEGORY = (
        "Unknown"  # placeholder for cells that do not have known labels aka STx cells
    )

    # Set the label column (SCANVI_CELLTYPE_KEY) to "Unknown" for all cells initially.
    # Ensures that STx cells are marked as unlabeled before we assign ref labels
    adata_combined.obs[SCANVI_CELLTYPE_KEY] = UNLABELED_CATEGORY

    logger.info(f"Columns in adata_combined: {list(adata_combined.obs.columns)}")

    # Ensure the batch column exists
    if BATCH_COL not in adata_combined.obs.columns:
        logger.error(f"Batch column '{BATCH_COL}' not found")
        return None, None

    # Assign real labels to reference cells
    ref_mask = adata_combined.obs[BATCH_COL] == REFERENCE_KEY

    # Ensure reference label column exists
    if REF_CELL_LABEL_COL not in adata_combined.obs.columns:
        logger.error(f"Cell type column '{REF_CELL_LABEL_COL}' not found")
        return None, None

    # Assign reference labels safely
    adata_combined.obs.loc[ref_mask, SCANVI_CELLTYPE_KEY] = adata_combined.obs.loc[
        ref_mask, REF_CELL_LABEL_COL
    ].astype(str)

    # Convert the column to categorical (recommended for scANVI)
    adata_combined.obs[SCANVI_CELLTYPE_KEY] = adata_combined.obs[
        SCANVI_CELLTYPE_KEY
    ].astype("category")

    logger.info(
        f"Reference cells with labels: {ref_mask.sum()}"
    )  # number of reference cells with labels.
    logger.info(
        f"Percent ref cells labeled: {100 * ref_mask.sum() / adata_combined.n_obs:.2f}%"
    )
    logger.info(
        f"Spatial cells (unlabeled): {(~ref_mask).sum()}"
    )  # number of spatial (unlabeled) cells.
    logger.info(
        f"Percent Spatial cells: {100 * (~ref_mask).sum() / adata_combined.n_obs:.2f}%"
    )
    logger.info(
        f"Unique cell types: {adata_combined.obs[SCANVI_CELLTYPE_KEY].value_counts()}"
    )

    logger.info("Initializing scANVI model...")
    SCANVI.setup_anndata(
        adata_combined,
        labels_key=SCANVI_CELLTYPE_KEY,
        unlabeled_category=UNLABELED_CATEGORY,
        batch_key=BATCH_COL,
    )
    scanvi_model = SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category=UNLABELED_CATEGORY,
        labels_key=SCANVI_CELLTYPE_KEY,
    )

    logger.info("Training scANVI model...")
    scanvi_model.train(
        max_epochs=MAX_EPOCHS_SCANVI,
        batch_size=1024,
        early_stopping=True,
        early_stopping_patience=15,
    )
    logger.info("scANVI training completed!")

    logger.info("Saving scANVI model...")
    scanvi_model.save(module_dir / "_scanvi_ref", overwrite=True)

    logger.info("Getting latent representation...")
    adata_combined.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation()

    logger.info("Visualizing scANVI latent space...")
    sc.pp.pca(adata_combined)
    sc.pp.neighbors(adata_combined, use_rep=SCANVI_LATENT_KEY)
    sc.tl.umap(adata_combined)
    sc.pl.umap(
        adata_combined,
        color=[BATCH_COL, REF_CELL_LABEL_COL],
        frameon=False,
        ncols=1,
        save=f"_{config.module_name}_scanvi_umap.png",
    )
    return adata_combined, scanvi_model


def extract_predictions_and_visualize(adata_combined, scanvi_model, adata, module_dir):
    """Extract predictions from scANVI model and copy them to original adata.

    Args:
    adata_combined : anndata.AnnData
        Combined reference and spatial data
    scanvi_model : scvi.model.SCANVI
        Trained scANVI model
    adata : anndata.AnnData
        Original spatial dataset to which predictions will be copied
    module_dir : Path
        Output directory for saving results

    Returns:
    adata : anndata.AnnData
        Original adata with scANVI predictions added
    """
    logger.info("Extracting scANVI predictions...")

    # Get predictions (hard labels)
    logger.info("Making predictions...")
    adata_combined.obs[SCANVI_LABEL_COL] = scanvi_model.predict(adata_combined)

    logger.info("Visualize predictions...")
    df = (
        adata_combined.obs.groupby([REF_CELL_LABEL_COL, SCANVI_LABEL_COL])
        .size()
        .unstack(fill_value=0)
    )
    norm_df = df / df.sum(axis=0)

    plt.figure(figsize=(8, 8))
    _ = plt.pcolor(norm_df)
    _ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
    _ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
    plt.xlabel("Predicted")
    plt.ylabel("Observed")
    plt.savefig(
        module_dir / "scanvi_confusion_matrix.png", dpi=300, bbox_inches="tight"
    )

    # Get prediction probabilities
    prediction_probs = scanvi_model.predict(adata_combined, soft=True)
    adata_combined.obsm["scanvi_probs"] = prediction_probs

    # Calculate prediction uncertainty (entropy)
    entropy = -np.sum(prediction_probs * np.log(prediction_probs + 1e-10), axis=1)
    adata_combined.obs["prediction_entropy"] = entropy

    # Calculate maximum probability (confidence)
    max_prob = np.max(prediction_probs, axis=1)
    adata_combined.obs["prediction_max_prob"] = max_prob

    # Extract spatial predictions
    spatial_mask = adata_combined.obs[BATCH_COL] == SPATIAL_KEY
    adata_spatial_predicted = adata_combined[spatial_mask].copy()

    logger.info(f"Spatial cells annotated: {adata_spatial_predicted.n_obs}")
    logger.info("Predicted cell type distribution:")
    logger.info(adata_spatial_predicted.obs[SCANVI_LABEL_COL].value_counts())

    # Map predictions back to original adata using cell names
    # Remove the suffix added by concat (index_unique="_")
    spatial_indices = adata_spatial_predicted.obs_names.str.replace(
        "_Spatial", "", regex=False
    )

    # Ensure indices match
    matching_mask = spatial_indices.isin(adata.obs_names)
    if not matching_mask.all():
        logger.warning(
            f"Some cell names don't match: {matching_mask.sum()}/{len(matching_mask)}"
        )

    # Create mapping
    cell_prediction_map = dict(
        zip(spatial_indices, adata_spatial_predicted.obs[SCANVI_LABEL_COL])
    )
    probability_map = dict(
        zip(spatial_indices, adata_spatial_predicted.obsm["scanvi_probs"])
    )
    entropy_map = dict(
        zip(spatial_indices, adata_spatial_predicted.obs["prediction_entropy"])
    )
    max_prob_map = dict(
        zip(spatial_indices, adata_spatial_predicted.obs["prediction_max_prob"])
    )

    # Copy predictions to original adata
    adata.obs[SCANVI_LABEL_COL] = adata.obs_names.map(cell_prediction_map)
    adata.obs["scANVI_prediction_probs"] = adata.obs_names.map(probability_map)
    adata.obs["scANVI_prediction_entropy"] = adata.obs_names.map(entropy_map)
    adata.obs["scANVI_prediction_confidence"] = adata.obs_names.map(max_prob_map)

    logger.info("scANVI predictions copied to original adata")

    # Create visualizations
    logger.info("Creating scANVI visualizations...")

    # UMAP for overview
    sc.pp.neighbors(adata_combined, use_rep=SCANVI_LATENT_KEY, n_neighbors=15)
    sc.tl.umap(adata_combined)

    # Plot results
    _, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Plot 1: UMAP by dataset
    sc.pl.umap(adata_combined, color=BATCH_COL, ax=axes[0, 0], show=False)
    axes[0, 0].set_title(BATCH_COL)

    # Plot 2: UMAP by cell type (reference) and predictions
    sc.pl.umap(adata_combined, color=SCANVI_LABEL_COL, ax=axes[0, 1], show=False)
    axes[0, 1].set_title("scANVI Predictions")

    # Plot 3: Prediction uncertainty
    sc.pl.umap(adata_combined, color="prediction_entropy", ax=axes[1, 0], show=False)
    axes[1, 0].set_title("Prediction Uncertainty (Entropy)")

    # Plot 4: Prediction confidence
    sc.pl.umap(adata_combined, color="prediction_max_prob", ax=axes[1, 1], show=False)
    axes[1, 1].set_title("Prediction Confidence (Max Probability)")

    plt.tight_layout()

    # Save plot
    plot_path = module_dir / "scanvi_integration_results.png"
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()

    logger.info(f"scANVI results plot saved: {plot_path}")

    return adata


def run_integration(
    config: IntegrateSCVIModuleConfig, io_config: IOConfig, base_config: Config
):
    """Integrate scRNAseq and STx data using scANVI and ingest.

    adata_ref_ingest and adata_ingest are used for ingest integration.
    adata_ref and data are used for scANVI integration.

    Args:
        config (IntegrateSCVIModuleConfig): Integration module configuration object.
        io_config (IOConfig): IO configuration object.
        base_config (Config): Base configuration object.

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

    # Set figure settings
    # configure_scanpy_figures(str(io_config.output_dir)) - not suitable fig settings

    # Log GPU availability
    if torch.cuda.is_available():
        logger.info(f"✓ GPU available: {torch.cuda.get_device_name(0)}")
        logger.info(f"✓ CUDA version: {torch.version.cuda}")
        logger.info("✓ scVI will use GPU acceleration")
    else:
        logger.info("No GPU detected - scVI will use CPU (slower training)")

    logger.info("Starting integration of scRNAseq and spatial transcriptomics data...")

    logger.info("Loading scRNAseq data from HLCA ...")
    adata_ref = sc.read_h5ad(ref_path)

    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(io_config.output_dir / "2_dimension_reduction" / "adata.h5ad")

    logger.info(
        "Checking if need to subsample reference dataset for faster integration..."
    )
    if adata_ref.n_obs > 500000:  # Only subsample very large references
        logger.info(f"Subsampling reference from {adata_ref.n_obs} to ~200k cells...")
        adata_ref_subset = subset_reference(
            base_config=base_config,
            adata_ref=adata_ref,
            target_cells=500000,
            stratify_col=REF_CELL_LABEL_COL,
            min_per_type=100,
            sampling_strategy="sqrt",
            rare_cell_boost=2.0,
        )
    else:
        logger.info(
            f"Reference dataset is small enough, using full dataset: {adata_ref.n_obs}"
        )
        adata_ref_subset = adata_ref

    # Clean up intermediate objects to free memory
    del adata_ref
    gc.collect()

    # Log available layers for debugging
    logger.info(f"Reference layers: {list(adata_ref_subset.layers.keys())}")
    logger.info(f"Spatial layers: {list(adata.layers.keys())}")

    # Log counts layer shapes (already verified to exist)
    logger.info(f"Reference counts shape: {adata_ref_subset.layers['counts'].shape}")
    logger.info(f"Spatial counts shape: {adata.layers['counts'].shape}")

    # Ensure REF_CELL_LABEL_COL column exists in spatial data
    if REF_CELL_LABEL_COL not in adata.obs.columns:
        adata.obs[REF_CELL_LABEL_COL] = "STx_UNKNOWN"  # placeholder for spatial cells

    logger.info("Verifying datasets for scVI integration...")
    scVI_integration_check(
        adata_ref_subset, batch_key=BATCH_COL, cell_type=REF_CELL_LABEL_COL
    )
    scVI_integration_check(adata, batch_key=BATCH_COL, cell_type=REF_CELL_LABEL_COL)

    logger.info("Preparing shared gene sets for scVI integration...")
    adata_ref_subset, adata_ingest = prepare_integrated_datasets(
        gene_id_dict_path, adata_ref_subset, adata
    )

    logger.info("Combine reference and query data ...")

    gc.collect()  # Force garbage collection before large operations

    # Create combined dataset with memory efficiency
    logger.info(
        f"Creating combined dataset from {adata_ref_subset.n_obs:,} ref + "
        f"{adata_ingest.n_obs:,} spatial cells"
    )
    adata_combined = anndata.concat(
        [adata_ref_subset, adata_ingest],
        join="inner",  # only keeps genes present in both datasets
        label=BATCH_COL,
        keys=[REFERENCE_KEY, SPATIAL_KEY],
        index_unique="_",
    )

    logger.info(
        f"Combined dataset created: {adata_combined.n_obs:,} cells x "
        f"{adata_combined.n_vars:,} genes"
    )

    # Verify counts layer was preserved during concat
    if "counts" not in adata_combined.layers:
        raise ValueError("Counts layer lost during concat operation!")
    else:
        logger.info("Counts layer preserved during concat")

    logger.info(f"Combined data layers: {list(adata_combined.layers.keys())}")

    # Clean up intermediate objects to free memory
    del adata_ref_subset
    del adata_ingest
    gc.collect()

    # Check dimensions
    logger.info(
        f"Combined dataset shape: {adata_combined.shape} "
        f"(cells x genes: {adata_combined.n_obs} x {adata_combined.n_vars})"
    )

    # Save combined data before integration
    combined_save_path = module_dir / "adata_combined_before_scVI.h5ad"
    adata_combined.write_h5ad(combined_save_path)
    logger.info(f"Combined data saved to {combined_save_path}")

    logger.info(
        "STEP 1: Harmonize scRNAseq reference dataset with STx dataset scVI model..."
    )
    adata_combined, trained_scvi_model = scVI_integration(
        config, adata_combined, module_dir
    )

    logger.info(f"scVI model training complete. adata combined: {adata_combined}")
    logger.info(f"scVI model training complete. Model: {trained_scvi_model}")

    logger.info("STEP 2. Transfer labels using scANVI model...")
    adata_combined, trained_scanvi_model = scANVI_label_transfer(
        config, adata_combined, trained_scvi_model, module_dir
    )

    # Extract scANVI predictions and copy to original adata
    logger.info("STEP 3: Extracting predicted labels from scANVI...")
    if adata_combined is not None and trained_scanvi_model is not None:
        logger.info("Visualize data following label transfer...")
        adata = extract_predictions_and_visualize(
            adata_combined, trained_scanvi_model, adata, module_dir
        )
    else:
        logger.error("scANVI integration failed. Skipping scANVI predictions.")

    logger.info("Saving integrated data...")
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Integrated data saved to {module_dir}.")

    logger.info("Integration complete.")
