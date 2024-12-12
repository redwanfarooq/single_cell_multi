# Description
Snakemake pipeline for merging, normalisation and integration of multisample/multimodal single cell datasets
- Modularised workflow can be modified and/or extended for different experiment designs
- Add as a submodule in a bioinformatics project GitHub repository
```
git submodule add https://github.com/redwanfarooq/single_cell_multi single_cell_multi
```
- Update submodule to the latest version
```
git submodule update --remote single_cell_multi
```

# Required software
1. Global environment
    - [Snakemake >=v7.31](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
    - [docopt >=v0.6](https://github.com/docopt/docopt)
    - [pandas >=v2.0](https://pandas.pydata.org/docs/getting_started/install.html)
    - [loguru >=v0.7](https://github.com/Delgan/loguru)
    - [h5py >=3.12](https://docs.h5py.org/en/latest/build.html)
2. Specific modules
    - [R >=v4.3](https://cran.r-project.org)
        * [docopt v0.7.1](https://CRAN.R-project.org/package=docopt)
        * [logger v0.3.0](https://CRAN.R-project.org/package=logger)
        * [qs v0.26.3](https://CRAN.R-project.org/package=qs)
        * [hdf5r v1.3.11](https://CRAN.R-project.org/package=hdf5r)
        * [tidyverse v2.0.0](https://CRAN.R-project.org/package=tidyverse)
        * [furrr v0.3.1](https://CRAN.R-project.org/package=furrr)
        * [MatrixExtra v0.1.15](https://CRAN.R-project.org/package=MatrixExtra)
        * [Seurat v5.1.0](https://CRAN.R-project.org/package=Seurat)
        * [Signac v1.14.0](https://CRAN.R-project.org/package=Signac)
        * [harmony v1.2.1](https://CRAN.R-project.org/package=harmony)
        * [Bioconductor v3.18](https://www.bioconductor.org/install/)
            + DropletUtils
            + batchelor
            + ensembldb
            + EnsDb.Hsapiens.v86
            + GenomicRanges
            + GenomeInfoDb
            + glmGamPoi
    - [MACS2 >=v2.2.9](https://github.com/macs3-project/MACS/wiki/Install-macs2)
    - [scvi-tools >=v1.2.0](https://scvi-tools.org/en/stable/install.html)

# Setup
1. Install software for global environment (requires Anaconda or Miniconda - see [installation instructions](https://conda.io/projects/conda/en/stable/user-guide/install/index.html))
    - Download [environment YAML](/resources/envs/snakemake.yaml)
    - Create new conda environment from YAML
    ```
    conda env create -f snakemake.yaml
    ```
2. Install software for specific module(s)
    - Manually install required software from source and check that executables are available in **PATH** (using `which`) *and/or*
    - Create new conda environments with required software from YAML (as above - download [environment YAMLs](/resources/envs)) *and/or*
    - Check that required software is available to load as environment modules (using `module avail`)
3. Set up pipeline configuration file **config/config.yaml** (see comments in file for detailed instructions)
4. Set up profile configuration file **profile/config.yaml** (see comments in file for detailed instructions)

# Run
1. Activate global environment
```
conda activate snakemake
```
2. Execute **run.py** in root directory

# Input
Pipeline requires the following input files/folders:

## General

**REQUIRED:**

1. Post-QC multimodal single cell count matrices in 10x-formatted HDF5 file for each sample
2. Cell barcode metadata table in delimited file format (e.g. TSV, CSV) for each sample
3. 10x-formatted indexed fragments file for each sample with ATAC data (if applicable)
4. Input metadata table in delimited file format (e.g. TSV, CSV) with the following required fields (with headers):
- **sample_id**: sample ID
- **hdf5**: path to 10x-formatted HDF5 file
- **metadata**: path to cell barcode metadata table
- **fragments**: path to 10x-formatted indexed fragments file (if applicable)
- **summits**: path to MACS2 peak summits BED file (if applicable)

# Output
Output directory will be created in specified locations with subfolders containing the output of each step specified in the module

# Modules

## Available modules
- default

## Adding new module
1. Add entry to module rule specifications file **config/modules.yaml** with module name and list of rule names
2. Add additional rule definition files in **modules/rules** folder (if needed)
- Rule definition file **must** also assign a list of pipeline target files generated by the rule to a variable with the same name as the rule
- Rule definition file **must** have the same file name as the rule with the file extension **.smk**
3. Execute **run.py** in root directory with `--update` flag (needs to be repeated if there are any further changes to the module rule specification in **config/modules.yaml**)