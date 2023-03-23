# Raw data for SERDP RDX project

## Amplicon sequence data

Raw sequencing data are not supplied in the GitHub repo and need to be supplied separately before running the pipeline.

This folder needs to contain the following six files comprising the SERDP RDX amplicon sequence dataset:

 - 16S_Undetermined_S0_L001_Index_TXRange_SERDP.fastq.gz (md5: ecc4f4bb44f1456248d08c25248e8123)
 - 16S_Undetermined_S0_L001_R1_TXRange_SERDP.fastq.gz (md5: ef6c4bb7868b2a84121325d185c2fd02)
 - 16S_Undetermined_S0_L001_R2_TXRange_SERDP.fastq.gz (md5: 269d9ee3de4cb07ea7df9d9295c196dd)
 - ITS_Undetermined_S0_L001_Index_TXRange_SERDP.fastq.gz (md5: 2668f2f8e874b4f416fb70c025679d9b)
 - ITS_Undetermined_S0_L001_R1_TXRange_SERDP.fastq.gz (md5: 3fa608d12006d323615c8e0d95d0e498)
 - ITS_Undetermined_S0_L001_R2_TXRange_SERDP.fastq.gz (md5: 6b6b6a032d83b69c5d790f4b5b3383c1)

Note that the pipeline takes gzipped files (*.fastq.gz) but these MD5 hashes are for the _uncompressed_ data. (i.e. use `zcat *.fastq.gz | md5sum` or the equivalent on your system to verify your checksums against the values given here.)

The files can be downloaded from: (add download instructions here once files are made available for public download)

The file names should match the raw_data files specified in [config.yaml]([/config/config.yaml).
