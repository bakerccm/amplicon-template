# Raw data for SERDP RDX project

See [scripts/download_data.sh](../scripts/download_data.sh) for shell code to download the sequence data.

## Amplicon sequence data

Raw sequencing data are not supplied in the GitHub repo and need to be supplied separately before running the pipeline.

This folder needs to contain the following files:

 - file1.fastq.gz (md5: md5md5md5md5md5md5md5md5md5md5md)
 - file2.fastq.gz (md5: md5md5md5md5md5md5md5md5md5md5md)
 - file3.fastq.gz (md5: md5md5md5md5md5md5md5md5md5md5md)
 - file4.fastq.gz (md5: md5md5md5md5md5md5md5md5md5md5md)
 - file5.fastq.gz (md5: md5md5md5md5md5md5md5md5md5md5md)
 - file6.fastq.gz (md5: md5md5md5md5md5md5md5md5md5md5md)

Note that the pipeline takes gzipped files (*.fastq.gz) but these MD5 hashes are for the _uncompressed_ data. (i.e. use `zcat *.fastq.gz | md5sum` or the equivalent on your system to verify your checksums against the values given here.)

The files can be downloaded from: (add download instructions here once files are made available for public download)

The file names should match the raw_data files specified in [config.yaml]([/config/config.yaml).
