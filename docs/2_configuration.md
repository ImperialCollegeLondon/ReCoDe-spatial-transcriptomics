# Configuration Guide

`recode_st` uses TOML configuration files to control which analysis modules to run and their parameters. Only modules defined in the configuration file will be executed.

## Basic Configuration Structure

```toml
# Global settings
log_level = "INFO"  # DEBUG, INFO, WARNING, ERROR, CRITICAL (default: INFO)
seed = 21122023     # Random seed for reproducibility

# IO config settings
[io]
base_dir = "."      # Base directory for input/output files (default: Current Working Directory)

# Module configurations (only include modules you want to run)
[modules.module_name]
module_name = "output_folder_name"
# module-specific parameters...
```

## IO Config

The IO config allows you to specify where to get data from and where to save it. The options
and their defaults are:

```toml
[io]
base_dir = "."                         # Base directory for the project (relative to the current working directory)
data_dir = "data"                      # Directory where the input data is stored (relative to base_dir)
output_dir = "analysis"                # Directory where the output will be saved (relative to base_dir)
xenium_dir = "xenium"                  # Directory containing Xenium data (relative to data_dir)
zarr_dir = "xenium.zarr"               # Directory where Zarr data will be stored (relative to data_dir)
area_path = "selected_cells_stats.csv" # Path to the CSV file containing selected cells statistics (relative to data_dir)
logging_path = "logs"                  # Directory where logs will be saved (relative to output_dir)
```

This resolves to the following default directory structure:

```bash
base_dir
├── analysis
│   ├── output_folder_name/
│   ├── ...
│   └── logs/
├── config.toml
└── data
    ├── selected_cells_stats.csv
    ├── xenium/
    └── xenium.zarr/
```

It is also possible to define absolute paths for each of these paths, which means you
are not locked into this particular directory structure.

## Available Modules

- `format_data` - Data formatting and preprocessing
- `quality_control` - Cell and gene filtering (requires `min_counts`, `min_cells`)
- `dimension_reduction` - PCA, UMAP, etc.
- `annotate` - Cell type annotation
- `view_images` - Tissue visualization (requires `gene_list`)
- `spatial_statistics` - Spatial analysis
- `muspan` - Advanced spatial analysis (requires MuSpAn license)
- `muspan_spatial_stat` - Spatial statistics and graph analysis (requires MuSpAn license)
- `muspan_spatial_graph` - Spatial graph construction (requires MuSpAn license)

## Example 1: Basic Quality Control Pipeline

```toml
log_level = "INFO"
seed = 12345

[io]
base_dir = "."

[modules.format_data]
module_name = "0_format"

[modules.quality_control]
module_name = "1_qc"
min_counts = 10
min_cells = 5

[modules.dimension_reduction]
module_name = "2_dr"
```

## Example 2: Full Analysis with Visualization and specifying the IO

```toml
log_level = "INFO"
seed = 12345

[io]
base_dir = "."                         # Base directory for the project (relative to the current working directory)
data_dir = "data"                      # Directory where the input data is stored (relative to base_dir)
output_dir = "analysis"                # Directory where the output will be saved (relative to base_dir)
xenium_dir = "xenium"                  # Directory containing Xenium data (relative to data_dir)
zarr_dir = "xenium.zarr"               # Directory where Zarr data will be stored (relative to data_dir)
logging_path = "logs"                  # Directory where logs will be saved (relative to output_dir)

[modules.format_data]
module_name = "0_format"

[modules.quality_control]
module_name = "1_qc"
min_counts = 15
min_cells = 3

[modules.dimension_reduction]
module_name = "2_dr"

[modules.annotate]
module_name = "3_annotate"

[modules.view_images]
module_name = "4_images"
gene_list = ["EPCAM", "CD3D", "CD68", "PTPRC"]

[modules.spatial_statistics]
module_name = "5_spatial"
```

## Running with Configuration

```bash
# Run with a configuration file
python -m recode_st config.toml
```

The pipeline will only run the modules you specify in your configuration file, allowing you to customize your analysis workflow.

## Developer Guide: Extending the Configuration System

This section explains how developers can add new configuration variables or create new modules in the ReCoDe Spatial Transcriptomics pipeline.

### Key Design Principles

- **All module configs inherit from `BaseModuleConfig`** which provides the required `module_name` field
- **Module configs use [Pydantic] models** for automatic validation and type checking
- **All modules are optional** in `ModulesConfig` (using `| None = None`)
- **Imports are lazy** - modules are only imported when they're actually used
- **Each module creates its own output directory** using `config.module_name`
- **Use type hints and docstrings** for all configuration parameters

This design ensures the pipeline remains modular, extensible, and efficient.

### Adding Variables to Existing Module Configs

To add new configurable parameters to an existing module:

1. **Update the module config class** in `src/recode_st/config.py`:

      ```python
      class QualityControlModuleConfig(BaseModuleConfig):
          """Configuration for the Quality Control module."""

          min_counts: int
          min_cells: int
          # Add your new parameter here
          new_parameter: float
          """Description of what this parameter does."""
      ```

1. **Update the module function** to use the new parameter:

      ```python
      def run_qc(config: QualityControlModuleConfig, io_config: IOConfig):
          # Use the new parameter
          threshold = config.new_parameter
          # ... rest of function
      ```

1. **Add the parameter to config files** like `config.toml`:

      ```toml
      [modules.quality_control]
      module_name = "1_qc"
      min_counts = 10
      min_cells = 5
      new_parameter = 0.5  # Add your new parameter
      ```

### Creating a New Module

To add a completely new analysis module:

1. **Create the module config class** in `src/recode_st/config.py`:

      ```python
      class MyNewModuleConfig(BaseModuleConfig):
          """Configuration for My New Module."""

          my_parameter: int
          """Description of this parameter."""

          another_parameter: tuple[str, ...]
          """List of items for analysis."""
      ```

1. **Add it to ModulesConfig** in the same file:

      ```python
      class ModulesConfig(BaseModel):
          # ...existing modules...

          my_new_module: MyNewModuleConfig | None = None
          """Configuration for My New Module."""
      ```

1. **Create the module file** `src/recode_st/my_new_module.py`:

      ```python
      """My new analysis module."""

      from logging import getLogger
      from recode_st.config import IOConfig, MyNewModuleConfig

      logger = getLogger(__name__)


      def run_my_new_module(config: MyNewModuleConfig, io_config: IOConfig):
          """Run my new analysis."""
          # Use config parameters
          param_value = config.my_parameter
          items = config.another_parameter

          # Create output directory
          module_dir = io_config.output_dir / config.module_name
          module_dir.mkdir(exist_ok=True)

          # Your analysis code here...
          logger.info(f"Running analysis with parameter: {param_value}")
      ```

1. **Add the conditional import** in `src/recode_st/__main__.py`:

      ```python
      def main(config: Config):
          # ...existing code...

          if config.modules.my_new_module:
              from recode_st.my_new_module import run_my_new_module

              logger.info("Running My New Module")
              run_my_new_module(config.modules.my_new_module, config.io)
      ```

1. **Add to configuration files**:

      ```toml
      [modules.my_new_module]
      module_name = "my_analysis"
      my_parameter = 42
      another_parameter = ["item1", "item2", "item3"]
      ```

[Pydantic]: https://docs.pydantic.dev/latest/
