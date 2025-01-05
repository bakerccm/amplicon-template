#!/bin/bash

# N.B. In principle this file should run when executed as a shell script.
# In practice individual commands often fail e.g. if the remote server
# rejects a request that comes too quickly after the last one, or if
# remote server is temporarily unavailable. It may be easier to run the
# commands individually.

# assumes starting from root of repository
    cd databases

# use with 16S and assignTaxonomy()
    wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
# use with 16S and assignTaxonomy()
    # same file as silva_nr99_v138.1_train_set.fa.gz above but going down to species level
    # wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz
# use with 16S and assignSpecies()
    wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz
# use with 16S and DECIPHER()
    wget http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData

# use with 16S and assignTaxonomy()
    wget https://zenodo.org/record/4310151/files/rdp_train_set_18.fa.gz
# use with 16S and assignSpecies()
    wget https://zenodo.org/record/4310151/files/rdp_species_assignment_18.fa.gz

# UNITE v10.0
# use with ITS and assignTaxonomy()
    wget https://s3.hpc.ut.ee/plutof-public/original/db5b0819-2766-4d2f-a2cb-9f100eb447cd.tgz

    # extract files from tarball
    tar -xf db5b0819-2766-4d2f-a2cb-9f100eb447cd.tgz

    # this extracts to two files: 
    #     sh_general_release_dynamic_s_04.04.2024_dev.fasta
    #     sh_general_release_dynamic_s_04.04.2024.fasta

    # recompress individual files
    gzip sh_general_release_dynamic_s_04.04.2024_dev.fasta # replaces with .fasta.gz
    gzip sh_general_release_dynamic_s_04.04.2024.fasta # replaces with .fasta.gz

