#!/bin/env python


"""
Generates wrapper script for Snakemake pipeline.
Requires:
    - Wrapper script template file
"""


# ==============================
# MODULES
# ==============================
import os
import string
import docopt
from loguru import logger


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC = """
Generate wrapper script for Snakemake pipeline

Usage:
  generate_wrapper.py --module=<module> --template=<template> [--file=<file>] [options]

Arguments:
  -m --module=<module>      Pipeline module (required)
  -t --template=<template>  Template file (required)
  -f --file=<file>          Output file name [default: Snakefile]

Options:
  -h --help                 Show this screen
"""


# ==============================
# GLOBAL VARIABLES
# ==============================
LOAD = """module {0}:
    snakefile: 'resources/modules/{0}.smk'
    config: config

use rule * from {0} as *"""


# ==============================
# FUNCTIONS
# ==============================
@logger.catch(reraise=True)
def _main(opt: dict) -> None:
    # Generate wrapper script
    logger.info("Generating pipeline wrapper script")
    logger.info("Template: {}", os.path.abspath(opt["--template"]))
    generate_wrapper(
        module=opt["--module"],
        template=opt["--template"],
        filename=opt["--file"],
    )
    logger.success("Output file: {}", os.path.abspath(opt["--file"]))


def generate_wrapper(module: str, template: str, filename: str | None = None) -> str:
    """
    Generate wrapper script for single cell data QC pipeline.

    Arguments:
        ``module``: Pipeline module.\n
        ``template``: Template file.\n
        ``filename``: Output file path or ``None``.

    Returns:
        Writes wrapper script to ``filename`` (if provided) and returns script as string.
    """
    with open(file=template, mode="r", encoding="UTF-8") as file:
        template = string.Template(file.read())

    out = template.substitute(LOAD=LOAD.format(module))

    if filename is not None:
        with open(file=filename, mode="w", encoding="UTF-8") as file:
            file.write(out)

    return out


# ==============================
# SCRIPT
# ==============================
if __name__ == "__main__":
    _main(opt=docopt.docopt(DOC))
