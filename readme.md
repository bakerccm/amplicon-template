# Code for processing and analyzing amplicon data

This repo contains code and metadata for analyzing amplicon sequencing data.

- The amplicon dataset is described [here](docs/readme_dataset.md).

- A workflow based on the [dada2](https://benjjneb.github.io/dada2/tutorial.html) pipeline is implemented in [Snakemake](https://snakemake.readthedocs.io) to process the raw sequence data. The workflow implementation is described in detail [here](docs/readme_pipeline.md).

## Reproducing our results

This pipeline runs on the ERDC Carpenter HPC cluster. It should be reproducible in other similar environments with relatively minimal modifications. Here are the basic steps you will need to follow.

### Ensure snakemake is installed

You need to have a suitable version of [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) installed. If you already have this from a previous project, there is no need to reinstall it. On some HPC systems, you may be able to skip installing it yourself and instead load an environment module that makes snakemake available to you.

If you do need to install snakemake, [conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) or (probably better) [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) are recommended for this. There is a yaml file at [envs/snakemake-7.25.0.yaml](/envs/snakemake-7.25.0.yaml) that should provide specifications for a suitable conda environment. Something like this should do the trick: 

```bash
mamba env create -f envs/snakemake-7.25.0.yaml
```

If you need to install conda/mamba, the Miniforge distribution suggested on the [mamba installation page](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) should work just fine.

### Clone the code repository

The code respository currently lives on the ERDC GitLab. Since you are reading this readme file, you have probably found it already. If you have not done so, use `git clone` to copy the repo to the cluster or other computer that you want to run it on.

### Add sequence data

The required raw sequence data files are described at [data/readme.md](/databases/readme.md). These files are too big to be stored in the code repository. (The sample [metadata](/metadata), however, are provided with the repo.) The sequence files need to be obtained separately and placed in the [data](/data) folder following the instructions in [data/readme.md](/data/readme.md). The checksums for the downloaded files should be compared against the values in [data/md5sums.txt](/data/md5sums.txt).

For now, the sequence data files are available in the Soil Micro shared drive at `datasets/...`. Eventually, these files will likely be made available for public download. Download details will be provided here, and shell commands to download the data at [scripts/download.sh](/scripts/download.sh).

### Add reference databases

Reference databases for taxonomy assignment should be located in [databases](/databases).

These databases are not included in the GitHub repo and should be downloaded according to the instructions in [databases/readme.md](/databases/readme.md). The checksums for the downloaded files should be compared against the values in [databases/md5sums.txt](/databases/md5sums.txt).

Shell commands to download the data are included at [scripts/download_databases.sh](/scripts/download_databases.sh).

### Install SEPP

Most of the required software can be downloaded and installed automatically by snakemake, but this is not the case for [SEPP](https://github.com/smirarab/sepp). Use the shell script at `scripts/install_sepp.sh` to install this software:

```bash
bash scripts/install_sepp.sh
```

### Use snakemake to create conda environments

The Carpenter HPC cluster has no internet connectivity from the compute nodes, so snakemake cannot create download software and create conda environments on the fly while running the pipeline. They will need to be installed first on the login nodes:

```bash
conda activate snakemake-7.25.0
snakemake -c 1 --use-conda --conda-create-envs-only all
```

You can skip this step if your computing environment permits internet connection during a compute job, or if you are supplying the software in some other way.

### Run pipeline

You should now be in a position to run the pipeline:

```bash
conda activate snakemake-7.25.0
snakemake -c 192 --use-conda all
```

If you are running this on a cluster, you probably want to run it in a non-interactive compute job. The script at [scripts/snakemake.pbs](/scripts/snakemake.pbs) provides an example for how this could be done on a PBS system.

You might also want to run part of the pipeline at a time in order to review results along the way. For example:

```bash
conda activate snakemake-7.25.0
snakemake -c 192 --use-conda all_demultiplex_qc
```

## Code reuse

The code in this repo is intended for general use with amplicon sequencing projects and is organized to facilitate reuse subject to the included [license](LICENSE). In particular, our implementation of the dada2 pipeline should be useable with similar 16S/ITS datasets by simply forking the repo and updating the parameters in the configuration file [config/config.yaml](/config/config.yaml) (e.g. names of input files, primer sequences, parameters used to control scripts, etc). However, the code may require customization depending on the specifics of the dataset -- e.g. file naming, file formatting, different genetic markers sequenced, input files organized differently, other processing steps such as decontamination, etc. The modular organization of the Snakemake workflow should facilitate such customization. If customization is required, the [snakefile](snakefile) and/or the [R code files](/code) should be updated as necessary.

---
Created by Chris Baker (@bakerccm)
