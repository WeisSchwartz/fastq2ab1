"""Logging and small helpers."""

from __future__ import annotations

import logging
import shutil
import sys


def configure_logging(verbose: bool = False) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
        force=True,
    )
    return logging.getLogger("fastq2ab1")


def require_binary(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise FileNotFoundError(
            f"Required external binary '{name}' was not found on PATH."
        )
    return path
