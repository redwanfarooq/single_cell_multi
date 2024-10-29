##########################################################################################
# Snakemake module
# Author: Redwan Farooq
# Module name: default
##########################################################################################


# --------------------------------------------------
# SETUP
# Load modules
import os
import yaml
from math import floor
from resources.scripts.rule import *

# Load and parse sample/library info from YAML file
with open(file=os.path.join(config.get("metadata_dir", "metadata"), "info.yaml"), mode="r", encoding="UTF-8") as file:
    info = yaml.load(stream=file, Loader=yaml.SafeLoader)
for key, value in parse_info(info).items():
    globals()[key] = value

# Set module rules list
module_rules = ['merge', 'normalise', 'integrate']

# Import rules
include: 'rules/merge.smk'
include: 'rules/normalise.smk'
include: 'rules/integrate.smk'

# Set targets list
targets = [x for rule in [merge, normalise, integrate] for x in rule]
# --------------------------------------------------


# --------------------------------------------------
# RULES
rule all:
	input: [os.path.abspath(x) for x in targets]
# --------------------------------------------------