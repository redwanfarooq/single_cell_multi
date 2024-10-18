#!/bin/bash -l
# Modified from default Snakemake jobscript.sh
# Using '#!/bin/bash -l' instead of '#!/bin/sh' to source /etc/profile at shell invocation (required to load environment modules)

# properties = {properties}
{exec_job}
