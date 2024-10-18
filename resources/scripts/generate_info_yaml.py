#!/bin/env python


"""
Generates info YAML for use with single cell multi pipeline.
Requires:
- Metadata table file with the following fields:
    sample_id: unique sample ID
    hdf5: path to 10x-formatted HDF5 file containing combined multimodal count matrices
    metadata: path to TSV files containing cell metadata
    fragments: path to to ATAC fragments file (if applicable)
"""


# ==============================
# MODULES
# ==============================
import os
import yaml
import math
import docopt
from loguru import logger
import pandas as pd


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC = """
Generate info YAML for use with single cell multi pipeline

Usage:
  generate_info_yaml.py --md=<md> --outdir=<outdir> [options]

Arguments:
  -m --md=<md>              Metadata table file (required)
  -o --outdir=<outdir>      Output directory (required)

Options:
  -h --help                 Show this screen
"""


# ==============================
# FUNCTIONS
# ==============================
@logger.catch(reraise=True)
def _main(opt: dict) -> None:
    # Read input CSV and check fields are valid
    md = pd.read_csv(opt["--md"], header=0, sep=None, engine="python")
    assert set(md.columns).issuperset(
        {"sample_id", "hdf5", "metadata"}
    ), "Invalid metadata table file."

    # Generate info YAML
    logger.info("Generating info YAML")
    generate_info_yaml(
        df=md,
        filename=os.path.join(opt["--outdir"], "info.yaml"),
    )
    logger.success(
        "Output file: {}", os.path.abspath(os.path.join(opt["--outdir"], "info.yaml"))
    )


def _filter_nan(d: dict) -> dict:
    return {
        k: _filter_nan(v) if isinstance(v, dict) else v
        for k, v in d.items()
        if not (isinstance(v, float) and math.isnan(v))
    }


def generate_info_yaml(df: pd.DataFrame, filename: str | None = None) -> dict:
    """
    Generate info YAML.

    Arguments:
        ``df``: DataFrame containing run metadata.\n
        ``filename``: Output file path or ``None``.

    Returns:
        Writes formatted YAML string to ``filename`` (if provided) and returns dictionary containing YAML data.
    """
    out = _filter_nan(df.set_index("sample_id").to_dict("index"))

    if filename is not None:
        with open(file=filename, mode="w", encoding="UTF-8") as file:
            yaml.dump(data=out, stream=file)

    return out


# ==============================
# SCRIPT
# ==============================
if __name__ == "__main__":
    _main(opt=docopt.docopt(DOC))
