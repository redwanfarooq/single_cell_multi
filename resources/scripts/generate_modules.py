#!/bin/env python


"""
Generates module scripts for Snakemake pipeline.
Requires:
    - Module rule specifications file
        YAML format with module names as keys and lists of module rules as values
    - Module script template file
"""


# ==============================
# MODULES
# ==============================
import os
import string
import yaml
import docopt
from loguru import logger


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC = """
Generate module scripts for Snakemake pipeline

Usage:
  generate_modules.py --modules=<modules> --template=<template> --outdir=<outdir> [options]

Arguments:
  -m --modules=<modules>    Module rule specifications file (required)
  -t --template=<template>  Template file (required)
  -o --outdir=<outdir>      Output directory (required)

Options:
  -h --help                 Show this screen
"""


# ==============================
# GLOBAL VARIABLES
# ==============================
LOAD = "include: 'rules/{}.smk'"
QUO = "'{}'"
UNQUO = "{}"


# ==============================
# FUNCTIONS
# ==============================
@logger.catch(reraise=True)
def _main(opt: dict) -> None:
    # Read input YAML
    with open(file=opt["--modules"], mode="r", encoding="UTF-8") as file:
        modules = yaml.load(stream=file, Loader=yaml.SafeLoader)

    # Generate module scripts
    logger.info("Updating module scripts")
    logger.info("Template: {}", os.path.abspath(opt["--template"]))
    for name, rules in modules.items():
        generate_module(
            name=name,
            rules=rules,
            template=opt["--template"],
            filename=os.path.join(opt["--outdir"], f"{name}.smk"),
        )
        logger.success(
            "Output file: {}",
            os.path.abspath(os.path.join(opt["--outdir"], f"{name}.smk")),
        )


def generate_module(
    name: str, rules: list[str], template: str, filename: str | None = None
) -> str:
    """
    Generate module script for preprocessing pipeline.

    Arguments:
        ``rules``: List of rules used in pipeline module.\n
        ``template``: Template file.\n
        ``filename``: Output file path or ``None``.

    Returns:
        Writes module script to ``filename`` (if provided) and returns script as string.
    """
    with open(file=template, mode="r", encoding="UTF-8") as file:
        template = string.Template(file.read())

    out = template.substitute(
        NAME=f"{name}",
        LOAD="\n".join([LOAD.format(x) for x in rules]),
        RULES_TO_STR=", ".join([QUO.format(x) for x in rules]),
        RULES_TO_VAR=", ".join([UNQUO.format(x) for x in rules]),
    )

    if filename is not None:
        with open(file=filename, mode="w", encoding="UTF-8") as file:
            file.write(out)

    return out


# ==============================
# SCRIPT
# ==============================
if __name__ == "__main__":
    _main(opt=docopt.docopt(DOC))
