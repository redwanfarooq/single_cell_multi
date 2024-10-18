#!/bin/env python


"""
Runs single cell multi pipeline.
"""


# ==============================
# MODULES
# ==============================
import os
import yaml
import docopt
from loguru import logger


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC = """
Run single cell multi pipeline

Usage:
  run.py [options]

Options:
  -h --help                 Show this screen
  -u --update               Update module scripts using rule specifications in 'config/modules.yaml' (will not run pipeline)
"""


# ==============================
# FUNCTIONS
# ==============================
@logger.catch(reraise=True)
def _main(opt: dict) -> None:
    # Get and execute shell command
    if not opt["--update"]:
        logger.info("Starting pipeline using module {}", MODULE)
    cmd = _get_cmd(update=opt["--update"])
    os.system(" && ".join(cmd))


def _cmd(*args):
    cmd = [" ".join(args)]
    return cmd


def _get_cmd(update: bool = False) -> list[str]:
    if update:
        cmd = _cmd(
            f"{SCRIPTS_DIR}/generate_modules.py",
            "--modules=config/modules.yaml",
            "--template=resources/templates/module.template",
            "--outdir=resources/modules",
        )
    else:
        cmd = _cmd(
            f"{SCRIPTS_DIR}/generate_wrapper.py",
            f"--module={MODULE}",
            "--template=resources/templates/wrapper.template",
        )
        cmd += _cmd(
            f"{SCRIPTS_DIR}/generate_info_yaml.py",
            f"--md={INPUT_TABLE}",
            f"--outdir={METADATA_DIR}",
        )
        cmd += _cmd("snakemake --profile=profile")
    return cmd


# ==============================
# SCRIPT
# ==============================
with open(file="config/config.yaml", mode="r", encoding="UTF-8") as file:
    config = yaml.load(stream=file, Loader=yaml.SafeLoader)
    SCRIPTS_DIR = config.get("scripts_dir", "resources/scripts")
    METADATA_DIR = config.get("metadata_dir", "metadata")
    try:
        INPUT_TABLE = os.path.join(METADATA_DIR, config["input"])
        MODULE = config["module"]
    except KeyError as err:
        raise KeyError(f"{err} not specified in '{file.name}'") from err


with open(file="config/modules.yaml", mode="r", encoding="UTF-8") as file:
    try:
        RULES = yaml.load(stream=file, Loader=yaml.SafeLoader)[MODULE]
    except KeyError as err:
        raise KeyError(f"Module {err} not specified in '{file.name}'") from err

if __name__ == "__main__":
    _main(opt=docopt.docopt(DOC))
