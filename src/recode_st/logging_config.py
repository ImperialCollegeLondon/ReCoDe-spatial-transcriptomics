"""Module for configuring logging in the recode_st package."""

import logging

from recode_st.paths import logging_path


def configure_logging(log_level=logging.INFO):
    """Configure logging for the package.

    Args:
        log_level: The logging level to set. Defaults to logging.INFO.
    """
    logger = logging.getLogger()
    logger.setLevel(log_level)

    ch = logging.StreamHandler()
    fh = logging.FileHandler(logging_path / "recode_st.log")

    formatter = logging.Formatter(
        "[%(asctime)s - %(name)s - %(levelname)s] %(message)s",
        datefmt="%Y/%m/%d %H:%M:%S",
    )

    ch.setFormatter(formatter)
    fh.setFormatter(formatter)

    logger.addHandler(ch)
    logger.addHandler(fh)
