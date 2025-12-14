"""The main entry point for the recode_st package."""

import sys
from logging import getLogger
from pathlib import Path

from recode_st.config import Config, load_config

logger = getLogger(__package__)


def main(config: Config):
    """Main function to run the recode_st package."""
    from recode_st.helper_function import configure_scanpy_figures, seed_everything
    from recode_st.logging_config import configure_logging

    configure_logging(config.io.logging_dir, config.log_level)

    # Apply Scanpy global visualization settings
    configure_scanpy_figures(config.io.output_dir)

    logger.info("Starting recode_st pipeline...")

    logger.info("Seeding everything...")
    seed_everything(config.seed)

    if config.modules.format_data:
        from recode_st.format_data import run_format

        logger.info("Running Module 0 - Format")
        run_format(config.io)

    if config.modules.subsample_data:
        from recode_st.subsample_data import run_subsampled_data

        logger.info("Running Optional Module - Subsample data")
        run_subsampled_data(config.modules.subsample_data, config.io, config)

    if config.modules.quality_control:
        from recode_st.qc import run_qc

        logger.info("Running Module 1 - Quality Control")
        run_qc(config.modules.quality_control, config.io)

    if config.modules.doublet_identification:
        from recode_st.doublet_identification import run_doublet_id

        logger.info("Running Module - Doublet identification")
        run_doublet_id(config.io, config.modules.doublet_identification)

    if config.modules.dimension_reduction:
        from recode_st.dimension_reduction import run_dimension_reduction

        logger.info("Running Module 2 - Dimension Reduction")
        run_dimension_reduction(config.modules.dimension_reduction, config.io)

    if config.modules.integrate_ingest:
        from recode_st.integrate_ingest import run_integration

        logger.info("Running Module - Integration: Ingest")
        run_integration(config.modules.integrate_ingest, config.io)

    if config.modules.integrate_scvi:
        from recode_st.integrate_scvi import run_integration

        logger.info("Running Module - Integration: scVI")
        run_integration(config.modules.integrate_scvi, config.io, config)

    if config.modules.annotate:
        from recode_st.annotate import run_annotate

        logger.info("Running Module - Annotate")
        run_annotate(config.modules.annotate, config.io)

    if config.modules.view_images:
        from recode_st.view_images import run_view_images

        logger.info("Running Module - View Images")
        run_view_images(config.modules.view_images, config.io)

    if config.modules.spatial_statistics:
        from recode_st.spatial_statistics import run_spatial_statistics

        logger.info("Running Module - Spatial Statistics")
        run_spatial_statistics(config.modules.spatial_statistics, config.io)

    if config.modules.drug2cell:
        from recode_st.drug2cell import run_drug2cell

        logger.info("Running Module - Drug2Cell")
        run_drug2cell(config.modules.drug2cell, config.io)

    if config.modules.muspan:
        from recode_st.muspan import run_muspan

        logger.info("Running Module - MuSpAn")
        run_muspan(config.modules.muspan, config.io)

    if config.modules.muspan_spatial_graph:
        from recode_st.ms_spatial_graph import run_muspan_graph

        logger.info("Running Module - MuSpAn Spatial Graph")
        run_muspan_graph(config.modules.muspan_spatial_graph, config.io)

    if config.modules.muspan_spatial_stat:
        from recode_st.ms_spatial_stat import run_muspan_stats

        logger.info("Running Module - MuSpAn Spatial Stats")
        run_muspan_stats(config.modules.muspan_spatial_stat, config.io)


if __name__ == "__main__":
    try:
        config_file = Path(sys.argv[1])
    except IndexError:
        raise Exception("No config file specified.")

    config = load_config(config_file)

    main(config)
