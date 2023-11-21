# Code for processing and analyzing amplicon data

This repo contains code and metadata for analyzing amplicon sequencing data.

- The amplicon dataset is described [here](docs/readme_dataset.md).

- A workflow based on the [dada2](https://benjjneb.github.io/dada2/tutorial.html) pipeline is implemented in [Snakemake](https://snakemake.readthedocs.io) to process the raw sequence data. The workflow implementation is described in detail [here](docs/readme_pipeline.md).

- Code for statistical analysis of the dataset following processing with dada2 is summarized [here](docs/readme_analysis.md).

## Reproducibility

To reproduce our analysis results, this GitHub repo needs to be augmented with our sequence data and the reference databases that we used.

### Sequence data

 Sequence data files should be located in [data](/data).

 The raw sequence data files are not included in the GitHub repo (although the sample [metadata](/metadata) are). The sequence data should be obtained following the instructions in [data/readme.md](/data/readme.md). The checksums for the downloaded files should be compared against the values in [data/md5sums.txt](/data/md5sums.txt).

 Shell commands to download the data are included at [scripts/download_data.sh](/scripts/download_data.sh).

### Reference databases

 Reference databases for taxonomy assignment should be located in [databases](/databases).

 These databases are not included in the GitHub repo and should be downloaded according to the instructions in [databases/readme.md](/databases/readme.md). The checksums for the downloaded files should be compared against the values in [databases/md5sums.txt](/databases/md5sums.txt).

 Shell commands to download the data are included at [scripts/download_databases.sh](/scripts/download_databases.sh).

## Code reuse

The code in this repo is intended for general use with amplicon sequencing projects and is organized to facilitate reuse subject to the included [license](LICENSE). In particular, our implementation of the dada2 pipeline should be useable with similar 16S/ITS datasets by simply forking the repo and updating the parameters in the configuration file [config/config.yaml](/config/config.yaml) (e.g. names of input files, primer sequences, parameters used to control scripts, etc). However, the code may require customization depending on the specifics of the dataset -- e.g. file naming, file formatting, different genetic markers sequenced, input files organized differently, other processing steps such as decontamination, etc. The modular organization of the Snakemake workflow should facilitate such customization. If customization is required, the [snakefile](snakefile) and/or the [R code files](/code) should be updated as necessary.

---
Maintained by [Chris Baker](https://github.com/bakerccm)
