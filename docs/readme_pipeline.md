# Amplicon data project

## Bioinformatics pipeline

### Overview

The pipeline in this repo is designed to process the six files (forward, reverse and index reads for the 16S and ITS datasets) provided by Argonne to CRREL prior to any demultiplexing or analysis. The output of the pipeline is a pair of `phyloseq` R objects (one for 16S and one for ITS) ready for statistical analysis, each comprising an ASV table, a table of sample metadata, and a table of taxonomic assignments.

This workflow is based on the [dada2 v1.18.0](https://benjjneb.github.io/dada2/tutorial.html) pipeline ([Callahan *et al*. 2016](#citations)).

The workflow is implemented here using the [Snakemake](https://snakemake.readthedocs.io) workflow engine ([Mölder *et al*. 2021](#citations)). This implementation provides a number of advantages including:

 - supports both 'dependencies' and 'reentrancy', *sensu* [Leipzig 2016](#citations);
 - supports flexible parallelization, which is especially useful for scaling the pipeline to larger datasets in a high-performance computing environment;
 - facilitates reproducibility by specifying software environments using [yaml files](/envs), and optionally automating the installation of those environments through Snakemake's interface with conda.

The pipeline runs with Snakemake v7.18.2 as specified in [envs/snakemake.yaml](/envs/snakemake.yaml).

### Input files

In addition to this repo, the pipeline requires raw sequence data and reference databases to be supplied in the folders [data](/data) and [databases](/databases).

The raw sequence data can be downloaded according to the instructions in [data/readme.md](/data/readme.md) and [data/download.sh](/data/download.sh), and checked against [data/md5sums.txt](/data/md5sums.txt). The databases used for the analysis can be downloaded according to the instructions in [databases/readme.md](/databases/readme.md) and [databases/download.sh](/databases/download.sh), and checked against [databases/md5sums.txt](/databases/md5sums.txt).

This pipeline may also be reused for other datasets. It assumes paired-end 16S and ITS sequencing data, with forward and reverse fastq files containing reads in matched order. Non-biological nucleotides (e.g. primers, adapters, linkers) should ideally be removed beforehand, although the `cutadapt` step of this workflow should hopefully clean up any residual primers. Metadata for the new dataset should be located at [metadata/samples_16S.tsv](/metadata/samples_16S.tsv) and [metadata/samples_ITS.tsv](/metadata/samples_ITS.tsv) in place of the project metadata.

### Configuration and computing environments

The configuration file at [config/config.yaml](/config/config.yaml) contains a variety of parameters that may be used to alter the running of the pipeline. The values in this file should be modified as required for other datasets (e.g. with input file names, primer sequences, and parameters to control execution of pipeline steps).

Software dependencies for the pipeline are recorded using the yaml files in [envs](/envs). These files may (optionally) be used by Snakemake to install any required software automatically using conda. Software may also be installed in other ways, including with the use of pre-installed environment modules if working in a high-performance computing environment. In this case, the contents of the yaml files may be used as a guide to the software versions that have been tested with this workflow.

### Pipeline steps

This section briefly outlines the steps in the workflow:

![rulegraph.png](/docs/rulegraph.png)

**barcode_list** extracts multiplexing barcodes from the metadata files supplied at [metadata/samples_16S.tsv](/metadata/samples_16S.tsv) and [metadata/samples_ITS.tsv](/metadata/samples_ITS.tsv) and saves them in the format required by the softward `iu-demultiplex`.

**demultiplex_16S** and **demultiplex_ITS** perform demultiplexing of the raw data using the software tool `iu-demultiplex`, creating one fastq file for each sample.

**compress_demultiplexed_fastq** compresses the demultiplexed output files using gzip to save storage, since `iu-demultiplex` does not output compressed files.

**fastqc** generates a quality report for each of the demultiplexed files using the software `fastqc`. It is not strictly required for the pipeline.

**filter_Ns** uses `R` to remove ambiguous bases from the demultiplexed fastq files, and checks the orientation of the primers supplied in [config/config.yaml](/config/config.yaml), in preparation for running `cutadapt`.

**cutadapt** uses the software `cutadapt` to find and trim primers from the demultiplexed fastq files. Searching is performed error-tolerantly using [default parameters](https://cutadapt.readthedocs.io), with the exception of the minimum sequence length. The latter is specified in the configuration file at [config/config.yaml](/config/config.yaml) and can be varied empirically to improve results. The minimum length should be >0 to avoid downstream problems in `dada`.

**read_quality_profiles** uses `R` to produce plots of sequence quality to help assess the output following `cutadapt`. It is not strictly required for the pipeline.

**filter_and_trim_16S** and **filter_and_trim_ITS** uses `R` to further filter and trim the output from `cutadapt` for the 16S and ITS datasets respectively. In both datasets, reads are truncated at the first instance of a quality score less than or equal to *truncQ*. After truncation, read1's and read2's are discarded if they have more than *maxEE_read1* or *maxEE_read2* expected errors, where expected errors are calculated from the quality scores as sum(10^(-Q/10)). In the 16S dataset, read1's and read2's are then truncated after *truncLen_read1* or *truncLen_read2* bases. Shorter reads are discarded entirely. This step is not performed for the ITS dataset, since length variation in ITS makes setting a hard cutoff inadvisable. Instead, a minimum length of *minLen* is imposed on the ITS dataset. All parameters are specified separately for each dataset in the configuration file at [config/config.yaml](/config/config.yaml) and can be varied empirically to improve results. The parameters *truncLen_read1* and *truncLen_read2* in particular should be chosen sensibly for the read/fragment length in the dataset (e.g. if this pipeline is employed on a dataset with longer reads, these values should not truncate reads too early).

**learnerrors_inputs** creates softlinks to the filtered and trimmed fastq files to group them by read direction in preparation for estimating error profiles.

**learnerrors** uses `dada2::learnErrors()` in `R` to estimate error profiles. The error profile is estimated using at least *nbases* bases, with this parameter configurable in [config/config.yaml](/config/config.yaml). Note that standard operation of learnErrors() involves reading in demultiplexed sequence files until at least *nbases* bases are read. Remaining sequence files are ignored. The bases are *not* drawn from across all files. One implication of this is that the result of this step depends on the order of the files supplied, though hopefully the error profile should be similar across all files in a run.

**dada** dereplicates each demultiplexed sequence file, uses the error profile calculated previously to infer the true composition of the sample using `dada2::dada()` in `R`, and finally merges paired ends. Note that this workflow implements this step in non-pooled mode, and may therefore discard some rare ASVs that would not be discarded if all samples were considered together. If rare taxa are of particular interest, it may be worth considering running this step in pooled or pseudo-pooled mode. Further details of these options are available [here](https://benjjneb.github.io/dada2/tutorial.html).

**dada_merge_samples** takes the output from **dada** for each of the demultiplexed samples and merges them into a single sequence table for each dataset.

**remove_chimeras** filters out chimeras using `dada2::removeBimeraDenovo()` in `R` using the default 'consensus' method.

**assign_taxonomy_16S** uses `dada2::assignTaxonomy()` to assign taxonomy given a database of reference sequences. For species-level assignments, two approaches are used: (i) only unambiguous identifications are allowed; and (ii) multiple species matches are allowed so long as they are consistent with the genus assigned by assignTaxonomy(). Results from both approaches are retained. Suitable reference databases should be supplied and specified in the configuration file at [config/config.yaml](/config/config.yaml). The RDP database and the Silva database are examples of suitable references, and details for obtaining them are provided in [databases/readme.md](/databases/readme.md).

**decipher** uses `DECIPHER` as an alternative method of assigning taxonomy to the 16S dataset.

**assign_taxonomy_ITS** uses `dada2::assignTaxonomy()` to assign taxonomy given a database of reference sequences. Suitable reference databases should be supplied and specified in the configuration file at [config/config.yaml](/config/config.yaml). The Unite database is an example of a suitable reference, and details for obtaining it is provided in [databases/readme.md](/databases/readme.md).

**export_to_phyloseq_16S** and **export_to_phyloseq_ITS** consolidate the processed sequence tables with read counts, the sample metadata, and the taxonomic assignments. These results are converted to `phyloseq` objects in `R` (one for each dataset). These objects are now ready for downstream statistical analysis using `phyloseq` or other software.

**remove_chloroplasts** applies several post-pipeline filters to the 16S dataset. It removes reads (i) designated as Eukaryotes; (ii) designated as Archaea; (iii) designated as mitochondria; (iv) shown as NA for Silva or RDP kingdom; (v) designated as chloroplasts. It also applies a length filter, only allowing through ASVs with sequences in the range 251:257 bp. The Silva, RDP and Decipher taxonomy results are jointly used for this filtering. See code in [remove_chloroplasts.R](/code/remove_chloroplasts.R) for details of how the filtering is executed.

**normalize_data** creates a copy of the processed `phyloseq` objects in which the read data are normalized using total sum scaling, i.e. dividing each read count by the total number of reads for the sample in the final dataset. For convenience, several versions of the normalized dataset are created in which results are aggregated at different taxonomic levels.

### Citations

Callahan, BJ, PJ McMurdie, MJ Rosen, AW Han, AJA Johnson and SP Holmes (2016). "DADA2: High-resolution sample inference from Illumina amplicon data." *Nature Methods* 13:581-583. doi:10.1038/nmeth.3869

Leipzig, Jeremy (2016). A review of bioinformatic pipeline frameworks. *Briefings in Bioinformatics* 2016:1–7. doi:10.1093/bib/bbw020

Mölder, F, KP Jablonski, B Letcher, MB Hall, CH Tomkins-Tinch, V Sochat, J Forster, S Lee, SO Twardziok, A Kanitz, A Wilm, M Holtgrewe, S Rahmann, S Nahnsen and J Köster (2021). Sustainable data analysis with Snakemake. *F1000Research* 10:33. doi:10.12688/f1000research.29032.1

---
Created by [Chris Baker](https://github.com/bakerccm)
