#!/bin/env python


"""
Runs single cell multi pipeline.
"""


# ==============================
# MODULES
# ==============================
import os
import glob
import sys
import shutil
import hashlib
from datetime import datetime
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
  -c --clean                Remove existing pipeline results in output directory
  -f --force                Overwrite existing pipeline results in output directory
  -u --update               Update module scripts using rule specifications in 'config/modules.yaml' (will not run pipeline)
"""


# ==============================
# FUNCTIONS
# ==============================
@logger.catch(reraise=True)
def _main(opt: dict) -> None:
    # Get and execute shell command
    cmd = _get_cmd(update=opt["--update"])
    if not opt["--update"]:
        # Check if output directory contains results of a previous pipeline run
        if os.path.exists(os.path.join(OUTPUT_DIR, ".pipeline", "md5sum")):
            if _file_to_str(
                os.path.join(OUTPUT_DIR, ".pipeline", "md5sum")
            ) != _get_hash(config, "VERSION", *METADATA):
                if opt["--clean"]:
                    logger.critical(
                        "Removing existing pipeline results in {}", OUTPUT_DIR
                    )
                    match str(
                        input("Are you sure you want to continue? (y/N) ")
                    ).lower():
                        case "y" | "yes":
                            shutil.rmtree(OUTPUT_DIR)
                        case _:
                            logger.error("Pipeline aborted")
                            sys.exit(1)
                elif opt["--force"]:
                    logger.critical(
                        "Overwriting existing pipeline results in {}", OUTPUT_DIR
                    )
                    match str(
                        input("Are you sure you want to continue? (y/N) ")
                    ).lower():
                        case "y" | "yes":
                            pass
                        case _:
                            logger.error("Pipeline aborted")
                            sys.exit(1)
                else:
                    logger.warning(
                        "Specified output directory contains results of a previous pipeline run from {} with different version, configuration and/or metadata ({}).\nUse --force to overwrite existing results (if required).\nUse --clean to remove existing results.",
                        _file_to_str(
                            os.path.join(OUTPUT_DIR, ".pipeline", "timestamp")
                        ),
                        os.path.join(OUTPUT_DIR, ".pipeline"),
                    )
                    sys.exit(0)
        logger.info("Starting pipeline using module {}", MODULE)
    if not os.system(" && ".join(cmd)) and not opt["--update"]:
        # Save run configuration and metadata to output directory for reproducibility
        os.makedirs(os.path.join(OUTPUT_DIR, ".pipeline"), exist_ok=True)
        with open(
            file=os.path.join(OUTPUT_DIR, ".pipeline", "md5sum"),
            mode="w",
            encoding="UTF-8",
        ) as f:
            f.write(_get_hash(config, "VERSION", *METADATA))
        with open(
            file=os.path.join(OUTPUT_DIR, ".pipeline", "timestamp"),
            mode="w",
            encoding="UTF-8",
        ) as f:
            f.write(str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        shutil.copy(src=CONFIG, dst=os.path.join(OUTPUT_DIR, ".pipeline"))
        shutil.copy(src="VERSION", dst=os.path.join(OUTPUT_DIR, ".pipeline"))
        for f in METADATA:
            if os.path.exists(f):
                shutil.copy(src=f, dst=os.path.join(OUTPUT_DIR, ".pipeline"))
        if os.path.exists("logs"):
            shutil.copytree(
                src="logs",
                dst=os.path.join(OUTPUT_DIR, ".pipeline", "logs"),
                dirs_exist_ok=True,
            )
        # Clean up cluster logs
        if glob.glob("sps-*"):
            os.system("rm -r sps-*")
        logger.success(
            "Pipeline completed successfully; run information saved to {}",
            os.path.join(OUTPUT_DIR, ".pipeline"),
        )


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


def _get_hash(options: dict, *args):
    x = [_file_to_str(_) if os.path.isfile(_) else "" for _ in args] if args else []
    x.append(yaml.dump(options, sort_keys=True))
    return hashlib.md5("".join(x).encode()).hexdigest()


def _file_to_str(path: os.PathLike) -> str:
    with open(file=path, mode="r", encoding="UTF-8") as f:
        return f.read().strip()


# ==============================
# SCRIPT
# ==============================
CONFIG = "config/config.yaml"
with open(file=CONFIG, mode="r", encoding="UTF-8") as file:
    config = yaml.load(stream=file, Loader=yaml.SafeLoader)
    SCRIPTS_DIR = config.get("scripts_dir", "resources/scripts")
    METADATA_DIR = config.get("metadata_dir", "metadata")
    try:
        INPUT_TABLE = os.path.join(METADATA_DIR, config["input"])
        OUTPUT_DIR = config["output_dir"]
        MODULE = config["module"]
    except KeyError as err:
        raise KeyError(f"{err} not specified in '{file.name}'") from err
METADATA = [_ for _ in [INPUT_TABLE] if _ is not None]


with open(file="config/modules.yaml", mode="r", encoding="UTF-8") as file:
    try:
        RULES = yaml.load(stream=file, Loader=yaml.SafeLoader)[MODULE]
    except KeyError as err:
        raise KeyError(f"Module {err} not specified in '{file.name}'") from err

if __name__ == "__main__":
    _main(opt=docopt.docopt(DOC))
