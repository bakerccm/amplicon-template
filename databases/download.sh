#!/bin/bash

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

# use with ITS and assignTaxonomy()
    wget https://files.plutof.ut.ee/public/orig/6A/F9/6AF94919CCB48307734D6256CACA50AE1ECBC0839F644D4B661E3673525E41A4.tgz
    tar -xf 6AF94919CCB48307734D6256CACA50AE1ECBC0839F644D4B661E3673525E41A4.tgz
    #
    # this extracts to a folder sh_general_release_s_10.05.2021, i.e.
    #
    # tree sh_general_release_s_10.05.2021
    # sh_general_release_s_10.05.2021
    # ├── sh_general_release_dynamic_s_10.05.2021_dev.fasta
    # └── sh_general_release_dynamic_s_10.05.2021.fasta
    #
    # recompress individual files
    cd sh_general_release_s_10.05.2021
    gzip sh_general_release_dynamic_s_10.05.2021_dev.fasta # replaces with .fasta.gz
    gzip sh_general_release_dynamic_s_10.05.2021.fasta # replaces with .fasta.gz
