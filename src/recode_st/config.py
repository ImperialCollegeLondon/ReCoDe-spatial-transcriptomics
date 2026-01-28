"""The configuration module for recode_st."""

from pathlib import Path
from typing import Literal, Self

import tomllib
from pydantic import BaseModel, DirectoryPath, Field, model_validator


class BaseModuleConfig(BaseModel):
    """Configuration for a specific module in recode_st."""

    module_name: str
    """The name of the module."""


class SubsampleModuleConfig(BaseModuleConfig):
    """Configuration for the Subsample module."""

    n_cells: int
    """Total number of cells to include in the subsampled dataset."""

    replace: bool = False
    """Whether to sample with replacement."""


class QualityControlModuleConfig(BaseModuleConfig):
    """Configuration for the Quality Control module."""

    subsample_data: bool = False
    """Whether to subsample the data for development purposes."""

    min_counts: int
    """Minimum number of counts required for a cell to pass filtering."""

    min_cells: int
    """Minimum number of cells expressed required for a gene to pass filtering"""

    min_cell_area: int
    """Minimum area to filter cells."""

    max_cell_area: int
    """Minimum area to filter cells."""

    remove_cells: list[str] = Field(default_factory=list)
    """List of cell types to remove from the dataset. Defaults to empty list."""

    norm_approach: Literal["scanpy_log", "sctransform", "cell_area", "none"] = (
        "cell_area"
    )
    """Normalization approach to use: 'scanpy_log' for Scanpy log normalization,
    'sctransform' for SCTransform normalization, or
    'cell_area' for normalization by cell area."""


class DoubletIdentificationModuleConfig(BaseModuleConfig):
    """Configuration for the Doublet Identification module."""


class DimensionReductionModuleConfig(BaseModuleConfig):
    """Configuration for the Quality Control module."""

    subsample_strategy: Literal["none", "compute", "load"] = "none"
    """Subsampling strategy for development dataset."""
    subsample_n_total: int = 5000
    """Total number of cells to include in the subsampled dataset."""
    subsample_min_cells_per_roi: int = 100
    """Minimum number of cells to include from each ROI in the subsampled dataset."""
    n_pca: int
    """number of principal components to compute"""
    n_neighbors: int
    """ number of neighbors for the neighborhood graph"""
    resolution: float
    """resolution for leiden clustering"""
    norm_approach: Literal["scanpy_log", "sctransform", "cell_area", "none"] = (
        "cell_area"
    )
    """Normalization approach to use."""
    obs_vis_list: list[str]
    """List of observation fields to visualize on UMAP."""
    marker_genes: list[str]
    """List of marker genes to visualize on UMAP."""


class IntegrateIngestModuleConfig(BaseModuleConfig):
    """Configuration for the Integration Ingest module."""


class IntegrateSCVIModuleConfig(BaseModuleConfig):
    """Configuration for the Integration scVI module."""


class AnnotateModuleConfig(BaseModuleConfig):
    """Configuration for the Annotate module."""

    cluster_name: str
    """Name of the cluster column in adata.obs to use for annotation."""

    clusters_label: str
    """Name of the new cluster column in adata.obs."""

    cluster_to_cell_type: dict[str, str]
    """Mapping from cluster labels to cell type annotations."""


class PsuedobulkModuleConfig(BaseModuleConfig):
    """Configuration for the Psuedobulk module."""

    annotation_var: str
    """Name of the annotation in adata.obs to use for pseudobulk aggregation."""

    subset_key: str
    """Key in adata.obs to subset the data for pseudobulk analysis."""


class ViewImagesModuleConfig(BaseModuleConfig):
    """Configuration for the View Images module."""

    gene_list: tuple[str, ...]
    """List of genes to visualize on tissue. A tuple of any length."""


class SpatialStatisticsModuleConfig(BaseModuleConfig):
    """Configuration for the Spatial Statistics module."""

    clusters_label: str
    """Name of the new cluster column in adata.obs."""


class Drug2CellModuleConfig(BaseModuleConfig):
    """Configuration for Drug2Cell."""

    drug_list: list[str]
    """List of drugs to visualize on UMAP and spatially in tissue"""


class MuspanModuleConfig(BaseModuleConfig):
    """Configuration for the Muspan module."""

    domain_name: str
    """Name of the domain for MuSpAn analysis."""

    transcripts_of_interest: tuple[str, ...]
    """List of transcripts of interest for MuSpAn analysis."""

    adata_cell_id: str
    """Column name in adata.obs for cell IDs."""

    cluster_labels: str
    """Column name in adata.obs for cell types."""


class MuspanSpatialStatModuleConfig(BaseModuleConfig):
    """Configuration for the Muspan Spatial Statistics module."""

    muspan_object: str
    """Name of the muspan object file."""

    cluster_labels: str
    """Name of the cluster labels field to use for analysis."""


class MuspanSpatialGraphModuleConfig(BaseModuleConfig):
    """Configuration for the Muspan Spatial Graph module."""

    muspan_object: str
    """Name of the muspan object file."""

    min_edge_distance: int
    """Minimum edge distance for spatial graphs."""

    max_edge_distance: int
    """Maximum edge distance for spatial graphs."""

    distance_list: tuple[int, ...]
    """List of distances for proximity networks."""

    min_edge_distance_shape: int
    """Minimum edge distance for shape-based networks."""

    max_edge_distance_shape: int
    """Maximum edge distance for shape-based networks."""

    k_list: tuple[int, ...]
    """List of k values for KNN networks."""


class FormatDataModuleConfig(BaseModuleConfig):
    """Configuration for the Format Data module."""


class ModulesConfig(BaseModel):
    """Configuration for all modules with optional fields."""

    format_data: FormatDataModuleConfig | None = None
    """Configuration for the Format Data module."""

    subsample_data: SubsampleModuleConfig | None = None
    """Configuration for the Subsample module."""

    quality_control: QualityControlModuleConfig | None = None
    """Configuration for the Quality Control module."""

    doublet_identification: DoubletIdentificationModuleConfig | None = None
    """Configuration for doublet detection with ovrlypy module."""

    dimension_reduction: DimensionReductionModuleConfig | None = None
    """Configuration for the Dimension Reduction module."""

    integrate_ingest: IntegrateIngestModuleConfig | None = None
    """Configuration for the Integration (ingest) module."""

    integrate_scvi: IntegrateSCVIModuleConfig | None = None
    """Configuration for the Integration (scVI) module."""

    annotate: AnnotateModuleConfig | None = None
    """Configuration for the Annotate module."""

    view_images: ViewImagesModuleConfig | None = None
    """Configuration for the View Images module."""

    spatial_statistics: SpatialStatisticsModuleConfig | None = None
    """Configuration for the Spatial Statistics module."""

    drug2cell: Drug2CellModuleConfig | None = None
    """Configuration for Drug2Cell module."""

    muspan: MuspanModuleConfig | None = None
    """Configuration for the Muspan module."""

    muspan_spatial_stat: MuspanSpatialStatModuleConfig | None = None
    """Configuration for the Muspan Spatial Statistics module."""

    muspan_spatial_graph: MuspanSpatialGraphModuleConfig | None = None
    """Configuration for the Muspan Spatial Graph module."""


class IOConfig(BaseModel):
    """Configuration for input and output directories."""

    base_dir: DirectoryPath = Path()
    """The base directory for all the input and output files."""

    data_dir: Path = Path("data")
    """The data directory for all the input files."""

    output_dir: Path = Path("out")
    """The output directory for all the output files."""

    xenium_dir: Path = Path("xenium")
    """The directory containing the Xenium input data."""

    output_data_dir: Path = Path("data")
    """The subdirectory under output_dir for processed data."""

    zarr_dir: Path = Path("xenium.zarr")
    """The directory containing the Zarr-formatted input data from each xenium ROI."""

    adata_dir: Path = Path("adata")
    """The directory containing the adata from each xenium ROI."""

    area_path: Path = Path("selected_cells_stats.csv")
    """The path to the CSV file containing selected cells statistics."""

    ref_path: Path  # no default, must be explicitly set in config
    """Path to the reference dataset in h5ad format."""

    gene_id_dict_path: Path  # no default, must be explicitly set in config
    """Path to the gene ID dictionary CSV file.
    The CSV should have 'ensembl_id' as the index column and a 'gene_symbol' column."""

    logging_dir: Path = Path("logs")
    """The directory for logging output."""

    @model_validator(mode="after")
    def resolve_paths(self) -> Self:
        """Resolve relative paths so they are absolute relative to base_dir."""
        # Input paths
        self.data_dir = self.base_dir / self.data_dir
        self.xenium_dir = self.data_dir / self.xenium_dir
        self.area_path = self.data_dir / self.area_path

        # Output paths
        self.output_dir = self.base_dir / self.output_dir
        self.output_data_dir = self.data_dir / self.output_data_dir
        self.zarr_dir = self.output_data_dir / self.zarr_dir
        self.adata_dir = self.output_data_dir / self.adata_dir
        self.logging_dir = self.output_dir / self.logging_dir

        # Ensure directories exists
        if not Path(self.data_dir).exists():
            self.data_dir.mkdir(parents=True, exist_ok=True)
        if not Path(self.output_dir).exists():
            self.output_dir.mkdir(parents=True, exist_ok=True)
        if not Path(self.logging_dir).exists():
            self.logging_dir.mkdir(parents=True, exist_ok=True)

        return self


class Config(BaseModel):
    """The possible configuration options for recode_st."""

    log_level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] = "INFO"
    """The logging level to use for the package."""

    seed: int
    """A random seed to use for reproducibility."""

    io: IOConfig = Field(default_factory=IOConfig)
    """Input and output directories configuration."""

    modules: ModulesConfig = Field(default_factory=ModulesConfig)
    """Configuration for all modules."""


def load_config(config_file: str | Path) -> Config:
    """Load the configuration from a TOML file.

    Args:
        config_file: The path to the configuration TOML file.

    Returns:
        An instance of Config with the loaded settings.
    """
    config_path = Path(config_file)
    if not config_path.exists():
        raise ValueError(f"Config file {config_path} not found.")

    with config_path.open("rb") as f:
        data = tomllib.load(f)

    return Config.model_validate(data)
