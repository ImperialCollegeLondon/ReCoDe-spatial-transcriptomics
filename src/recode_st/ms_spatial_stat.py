"""Muspan module - spatial statistics and graph analysis."""

import warnings
from logging import getLogger
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from recode_st.config import IOConfig, MuspanSpatialStatModuleConfig
from recode_st.helper_function import configure_scanpy_figures

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)


def calculate_and_plot_cross_pcf(
    ms,
    domain,
    cluster_labels,
    module_dir,
    unique_clusters,
    cell_type_indices=(0, 1),
    max_R=200,
    annulus_step=5,
    annulus_width=25,
    visualise_output=True,
):
    """Calculates and plots the cross-PCF for two selected cell types."""
    cell_type_indices = list(cell_type_indices)

    try:
        cell_type_1 = unique_clusters[cell_type_indices[0]]
        cell_type_2 = unique_clusters[cell_type_indices[1]]
    except IndexError:
        raise ValueError(
            "Cell type indices are out of range of the unique cluster list."
        )

    logger.info(f"Calculating cross-PCF for {cell_type_1} and {cell_type_2}...")

    # Query populations
    pop_A = ms.query.query(domain, ("label", cluster_labels), "is", str(cell_type_1))
    pop_B = ms.query.query(domain, ("label", cluster_labels), "is", str(cell_type_2))

    # Calculate cross-PCF
    r, PCF = ms.spatial_statistics.cross_pair_correlation_function(
        domain=domain,
        population_A=pop_A,
        population_B=pop_B,
        max_R=max_R,
        annulus_step=annulus_step,
        annulus_width=annulus_width,
        visualise_output=visualise_output,
    )

    # Save PCF plot
    pcf_plot_path = (
        module_dir / f"cross_pair_correlation_function_{cell_type_1}_{cell_type_2}.png"
    )
    plt.savefig(pcf_plot_path)
    logger.info(f"Cross-PCF plot saved at {pcf_plot_path}")

    # Visualize and save cell type points
    query_1_2 = ms.query.query(
        domain,
        ("label", cluster_labels),
        "in",
        [str(cell_type_1), str(cell_type_2)],
    )

    fig, ax = ms.visualise.visualise(
        domain, color_by=("label", cluster_labels), objects_to_plot=query_1_2
    )

    ax.set_title(f"{cell_type_1} vs {cell_type_2}", fontsize=15)
    ax.tick_params(axis="both", which="major", labelsize=10)
    ax.set_xlabel(f"{cell_type_1}", fontsize=15)
    ax.set_ylabel(f"{cell_type_2}", fontsize=15)

    vis_plot_path = module_dir / f"visualize_{cell_type_1}_{cell_type_2}.png"
    plt.savefig(vis_plot_path)
    logger.info(f"Visualization saved at {vis_plot_path}")


def calculate_pairwise_cross_pcf(
    ms,
    domain,
    module_dir,
    cluster_labels: str,
    unique_clusters: list[str],
):
    """Calculates and plots the cross-PCF for all unique unordered cell type pairs."""
    num_clusters = len(unique_clusters)
    fig, axes = plt.subplots(num_clusters, num_clusters, figsize=(40, 40))

    for i, cluster_i in enumerate(unique_clusters):
        for j, cluster_j in enumerate(unique_clusters):
            if j <= i:
                axes[i, j].axis("off")  # Optional: remove lower triangle plots
                continue

            pop_A = ms.query.query(domain, ("label", cluster_labels), "is", cluster_i)
            pop_B = ms.query.query(domain, ("label", cluster_labels), "is", cluster_j)
            logger.info(f"Calculating cross-PCF: {cluster_i} vs {cluster_j}")

            r, pcf = ms.spatial_statistics.cross_pair_correlation_function(
                domain,
                pop_A,
                pop_B,
                max_R=200,
                annulus_step=5,
                annulus_width=25,
            )

            ax = axes[i, j]
            ax.plot(r, pcf)
            ax.axhline(2, color="k", linestyle=":")
            ax.tick_params(axis="both", which="major", labelsize=15)
            ax.set_ylabel(f"$g_{{{cluster_i},{cluster_j}}}(r)$", fontsize=20)
            ax.set_xlabel("$r$", fontsize=20)

    plt.tight_layout()
    output_path = Path(module_dir) / "cross_pair_correlation_function_all.png"
    plt.savefig(output_path)
    logger.info(f"Cross-PCF matrix plot saved at {output_path}")


def run_muspan_stats(config: MuspanSpatialStatModuleConfig, io_config: IOConfig):
    """Run Muspan spatial statistics analysis on Xenium data."""
    try:
        import muspan as ms
    except ModuleNotFoundError as err:
        logger.error(
            "Could not load MuSpAn. Install with:\n"
            "    pip install 'recode_st[muspan]' @ git+https://github.com/ImperialCollegeLondon/ReCoDe-spatial-transcriptomics.git"
        )
        raise err

    module_dir = io_config.output_dir / config.module_name
    module_dir.mkdir(exist_ok=True)

    # Set figure settings to ensure consistency across all modules
    configure_scanpy_figures(str(io_config.output_dir))
    cmap = sns.color_palette("Spectral", as_cmap=True)

    # Load MuSpAn object
    logger.info("Loading MuSpAn object...")
    domain = ms.io.load_domain(path_to_domain=str(module_dir / config.muspan_object))

    # Get cluster labels
    cluster_labels = config.cluster_labels
    all_cluster_labels = domain.labels[cluster_labels]["labels"].tolist()
    unique_clusters = np.unique(all_cluster_labels).astype(str).tolist()
    logger.info(f"Found {len(unique_clusters)} unique cell types.")

    # Run pairwise analysis for specific pair
    logger.info("Calculating pairwise cross-PCF for selected cell types...")
    calculate_and_plot_cross_pcf(
        ms=ms,
        domain=domain,
        cluster_labels=cluster_labels,
        module_dir=module_dir,
        unique_clusters=unique_clusters,
        cell_type_indices=(0, 1),
    )

    # Full pairwise matrix
    logger.info("Calculating pairwise cross-PCF for all cell types...")
    calculate_pairwise_cross_pcf(
        ms=ms,
        domain=domain,
        module_dir=module_dir,
        cluster_labels=cluster_labels,
        unique_clusters=unique_clusters,
    )
