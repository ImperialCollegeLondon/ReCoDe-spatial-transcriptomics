"""Module containing helper functions."""

import os
import random

import numpy as np
import scanpy as sc
import seaborn as sns
import torch


def seed_everything(seed: int):
    """Set random seed on every random module for reproducibility.

    Args:
        seed: The seed value to set for random number generation.
    """
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = True
    elif torch.backends.mps.is_available():
        torch.mps.manual_seed(seed)
    else:
        pass


def configure_scanpy_figures(figdir: str, dpi: int = 300, fontsize: int = 16):
    """Configure global Scanpy visualization settings.

    This helper centralizes all figure-related parameters so that Scanpy
    produces consistent, publication-quality plots. It also sets the directory
    where Scanpy saves figures and returns commonly used visualization assets
    (e.g., colormaps).

    Args:
        figdir (str):
            Directory where Scanpy should save all generated figures.
        dpi (int, optional):
            Resolution to use for both on-screen display and saved figures.
            Defaults to 300.
        fontsize (int, optional):
            Base font size for all Scanpy plots. Defaults to 16.

    Returns:
        dict: A dictionary containing visualization resources. Currently
        includes:
            - ``"cmap"``: A seaborn colormap suitable for continuous variables.
    """
    sc.set_figure_params(
        dpi=dpi,
        dpi_save=dpi,
        frameon=False,
        vector_friendly=True,
        fontsize=fontsize,
        facecolor="white",
        figsize=(6, 4),
    )

    cmap = sns.color_palette("crest", as_cmap=True)

    return {"cmap": cmap}
