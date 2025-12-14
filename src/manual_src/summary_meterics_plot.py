"""Script to read and visualize the combined metrics summary CSV file."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Setup
base_dir = Path(
    "/Volumes/sep22/home/wet_lab/_Experiments/009_ST_Xenium/data/out_data/summary_metrics/"
)

out_dir = base_dir / "figs"
out_dir.mkdir(parents=True, exist_ok=True)

# Load combined dataframe
df = pd.read_excel(base_dir / "csv/combined_run_dfs.xlsx")
print(df.columns)

# Color palettes
condition_color_dict = {"PM08": "#BC3C29", "IPF": "#E18727", "COPD": "#0072B5"}
treatment_color_dict = {"Sham": "#FFDC91", "Treamtment": "#20854E"}

# Global style
sns.set_theme(style="white")  # removes grid lines
plt.rcParams.update(
    {
        "figure.figsize": (6, 6),
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "axes.titlepad": 12,
        "axes.labelpad": 8,
        "savefig.bbox": "tight",
        "savefig.dpi": 300,
        "figure.autolayout": True,
    }
)

# Data prep
df["Condition"] = pd.Categorical(
    df["Condition"], categories=["PM08", "COPD", "IPF"], ordered=True
)
df["norm_transcript"] = (
    df["total_high_quality_decoded_transcripts"] / df["total_cell_area"]
)
df["norm_cells"] = df["num_cells_detected"] / df["total_cell_area"]


# Helper function
def plot_violin(data, x, y, palette, title, ylabel, filename, order=None, ylim=None):
    """Plot violin plots.

    Args:
        data (_type_): _description_
        x (_type_): _description_
        y (_type_): _description_
        palette (_type_): _description_
        title (_type_): _description_
        ylabel (_type_): _description_
        filename (_type_): _description_
        order (_type_, optional): _description_. Defaults to None.
        ylim (_type_, optional): _description_. Defaults to None.
    """
    plt.figure()
    sns.violinplot(
        x=x,
        y=y,
        data=data,
        inner="point",
        palette=palette,
        order=order,
    )
    plt.title(title)
    plt.xlabel(x)
    plt.ylabel(ylabel)
    if ylim:
        plt.ylim(ylim)
    plt.savefig(filename)
    plt.close()


def raincloud_plot(
    data,
    x,
    y,
    hue,
    order=None,
    hue_order=None,
    palette=None,
    box_width=0.15,
    violin_width=0.6,
    jitter=0.15,
    figsize=(8, 5),
    title=None,
    ylabel=None,
    xlabel=None,
    ylim=None,
    filename=None,
):
    """Create a raincloud plot.

    Half violin (distribution), overlaid boxplot (summary), and jittered points
    colored by `hue` to show subgroup structure (e.g. run effects).

    Args:
        data: DataFrame containing the variables to plot.
        x: Categorical variable for the x-axis (e.g., Condition).
        y: Numeric variable to visualize.
        hue: Subgroup/category used to color individual points (e.g., run).
        order: Optional sequence defining the order of x-axis categories.
        hue_order: Optional sequence defining the order of hue categories.
        palette: Color mapping for hue categories.
        box_width: Width of the central boxplot element.
        violin_width: Width of the half-violin distribution shape.
        jitter: Horizontal jitter magnitude for point spreading.
        figsize: Figure size in inches (width, height).
        title: Optional plot title.
        ylabel: Optional y-axis label.
        xlabel: Optional x-axis label.
        ylim: Optional y-axis limits (tuple).
        filename: Optional path to save the figure. If None, plot is not saved.

    Returns:
        None
    """
    plt.figure(figsize=figsize)

    # ----------------------------------------------------------------------
    # 1. HALF VIOLIN (left half)
    # We use a full violin but clip to one side using a custom helper.
    # ----------------------------------------------------------------------
    sns.violinplot(
        data=data,
        x=x,
        y=y,
        order=order,
        inner=None,
        linewidth=0,
        palette=None,  # no coloring here; just outline disabled
        cut=0,
        width=violin_width,
        color="#D1EEEE",
    )

    # Clip the violins to left side
    ax = plt.gca()
    for collection in ax.collections:
        paths = collection.get_paths()
        for p in paths:
            v = p.vertices
            v[:, 0] = np.minimum(v[:, 0], np.median(v[:, 0]))  # keep the left half

    # ----------------------------------------------------------------------
    # 2. BOX PLOTS (centered)
    # ----------------------------------------------------------------------
    sns.boxplot(
        data=data,
        x=x,
        y=y,
        order=order,
        width=box_width,
        showcaps=True,
        boxprops={"facecolor": (1, 1, 1, 0.5), "zorder": 10},
        showfliers=False,
        whiskerprops={"linewidth": 1},
        medianprops={"color": "black", "linewidth": 2},
    )

    # ----------------------------------------------------------------------
    # 3. JITTERED POINTS (colored by run)
    # ----------------------------------------------------------------------
    sns.stripplot(
        data=data,
        x=x,
        y=y,
        hue=hue,
        order=order,
        hue_order=hue_order,
        palette=palette,
        dodge=True,
        jitter=jitter,
        alpha=0.8,
        zorder=20,
    )

    # ----------------------------------------------------------------------
    # 4. Labels and formatting
    # ----------------------------------------------------------------------
    if xlabel:
        plt.xlabel(xlabel)
    else:
        plt.xlabel(x)

    if ylabel:
        plt.ylabel(ylabel)
    else:
        plt.ylabel(y)

    if title:
        plt.title(title)

    if ylim:
        plt.ylim(ylim)

    plt.legend(title=hue, bbox_to_anchor=(1.02, 1), loc="upper left")

    plt.tight_layout()

    if filename:
        plt.savefig(filename, dpi=300)
        plt.close()


# Plots
plots = [
    (
        "total_cell_area",
        "Total Cell Area by Condition",
        "Total Cell Area (um²)",
        None,
        "total_cell_area_by_condition.png",
    ),
    (
        "norm_transcript",
        "High Quality Decoded Transcripts Normalized by Cell Area (um²)",
        "Normalized Transcripts (Mean Q >20)",
        None,
        "norm_high_Q_transcript_condition.png",
    ),
    (
        "fraction_transcripts_decoded_q20",
        "Fraction of Transcripts Decoded (Mean Q >20) by Condition",
        "Fraction of Transcripts Decoded (Mean Q >20)",
        (0, 1),
        "fraction_transcripts_decoded_q20_condition.png",
    ),
    (
        "norm_cells",
        "Number of Cells Detected Normalized by Cell Area (um²) by Condition",
        "Normalized Number of Cells Detected (Cells/um²)",
        None,
        "norm_cells_detected_condition.png",
    ),
    (
        "fraction_empty_cells",
        "Fraction of Empty Cells by Condition",
        "Fraction of Empty Cells",
        (-0.001, 0.5),
        "fraction_empty_cells_condition.png",
    ),
    (
        "median_genes_per_cell",
        "Median Genes per Cell by Condition",
        "Median Genes per Cell",
        None,
        "median_genes_per_cell_condition.png",
    ),
    (
        "median_transcripts_per_cell",
        "Median Transcripts per Cell by Condition",
        "Median Transcripts per Cell",
        None,
        "median_transcripts_per_cell_condition.png",
    ),
]

for y, title, ylabel, ylim, filename in plots:
    plot_violin(
        data=df,
        x="Condition",
        y=y,
        palette=condition_color_dict,
        title=title,
        ylabel=ylabel,
        filename=out_dir / filename,
        order=["PM08", "COPD", "IPF"],
        ylim=ylim,
    )

    raincloud_plot(
        data=df,
        x="Condition",
        y=y,
        hue="run",
        order=["PM08", "COPD", "IPF"],
        hue_order=sorted(df["run"].unique()),
        palette="Set2",
        title=title,
        ylabel=ylabel,
        xlabel="Condition",
        ylim=ylim,
        filename=out_dir / f"raincloud_{filename}",
    )

print("Plots saved successfully.")
