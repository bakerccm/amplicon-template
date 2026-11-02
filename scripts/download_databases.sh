#!/bin/bash

# assumes starting from root of repository
    cd databases

# use with 16S and assignTaxonomy()
    wget https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz
# use with 16S and assignTaxonomy()
    # same file as silva_nr99_v138.1_train_set.fa.gz above but going down to species level
    # wget https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toSpecies_trainset.fa.gz
# use with 16S and assignSpecies()
    wget https://zenodo.org/records/14169026/files/silva_v138.2_assignSpecies.fa.gz
    
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
