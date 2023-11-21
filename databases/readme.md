# Taxonomic Databases

See [scripts/download_databases.sh](../scripts/download_databases.sh) for code to download databases.

General information on training sets for use with dada2: https://benjjneb.github.io/dada2/training.html

## Silva 138.1 prokaryotic SSU taxonomic training data formatted for DADA2

See [this page](https://zenodo.org/record/4587955#.YMEPdC1r3UI) for notes and download links.

  - Use with 16S & assignTaxonomy(): [silva_nr99_v138.1_train_set.fa.gz](https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz) (md5:6b41db7139834c71171f8ce5b5918fc6)

  The alternative [silva_nr99_v138.1_wSpecies_train_set.fa.gz](https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz) (md5:ba13ab369161e0cb85df7e0ee3a4182e) is the same sequences but with taxonomy going down to species. But we are assigning species-level taxonomy with addSpecies()/assignSpecies() instead.

  - Use with 16S & assignSpecies(): [silva_species_assignment_v138.1.fa.gz](https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz) (md5:f21c2d97c79ff07c17949a9622371a4c)

  - Use with 16S & DECIPHER (modified SILVA SSU r138 dataset): [SILVA_SSU_r138_2019.RData](http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData) (md5:cb983b6a5e8cdb46f8c88b5afae21f66)

  Citations:

  Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucl. Acids Res. 41 (D1): D590-D596. [DOI:10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)

  Yilmaz P, Parfrey LW, Yarza P, Gerken J, Pruesse E, Quast C, Schweer T, Peplies J, Ludwig W, Glöckner FO (2014) The SILVA and "All-species Living Tree Project (LTP)" taxonomic frameworks. Nucl. Acids Res. 42:D643-D648 [DOI:10.1093/nar/gkt1209](https://doi.org/10.1093/nar/gkt1209)

  Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. 2016. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13:581–583. [DOI:10.1038/nmeth.3869](http://doi.org/10.1038/nmeth.3869)
  
  McLaren, Michael R., & Callahan, Benjamin J. (2021). Silva 138.1 prokaryotic SSU taxonomic training data formatted for DADA2 [Data set]. Zenodo. [DOI:10.5281/zenodo.4587955](http://doi.org/10.5281/zenodo.4587955)

## RDP taxonomic training data formatted for DADA2 (RDP trainset 18/release 11.5)

See [this page](https://zenodo.org/record/4310151#.YMEX9C1r3UI) for notes and download links.

  - Use with 16S & assignTaxonomy(): [rdp_train_set_18.fa.gz](https://zenodo.org/record/4310151/files/rdp_train_set_18.fa.gz) (md5:693c1ade0caf223e1fb7c6cf4da6f272)

  - Use with 16S & assignSpecies(): [rdp_species_assignment_18.fa.gz](https://zenodo.org/record/4310151/files/rdp_species_assignment_18.fa.gz) (md5:1a465122bdf670a74eb8873de042bad6)
  
  Citation:
  
  Callahan, Benjamin. (2020). RDP taxonomic training data formatted for DADA2 (RDP trainset 18/release 11.5) [Data set]. Zenodo. [DOI:10.5281/zenodo.4310151](http://doi.org/10.5281/zenodo.4310151)
  
## UNITE ITS database (Version 9.0 2023-07-25)

See [this page](https://unite.ut.ee/repository.php) for notes and download links (use general FASTA releases with DADA2)

  Fungi, includes global and 97% singletons: [https://doi.org/10.15156/BIO/2938068](https://doi.org/10.15156/BIO/2938068)

  - Use with ITS and assignTaxonomy(): [download](https://files.plutof.ut.ee/public/orig/F6/7D/F67DDF951A780D041DE915EBC6A0C006C1F8E9ECB7F3158BC8FEB32963CFF748.tgz) (md5: 4ceb511e0376295b7dad55b7e723367b) This downloads as a tar archive (.tgz) which extracts to two files: sh_general_release_dynamic_s_25.07.2023_dev.fasta (md5: 2e1789cad6e15a88d518a06f60ad1b47) and sh_general_release_dynamic_s_25.07.2023.fasta (md5: 3c7f6b7e826a3fe1de5ad06f128fc556). The 'dev' version apparently contains sequence flanking ITS which is removed in the other fasta file. Some classifiers (e.g. the old QIIME2 classifier) are better trained on the full length sequences. May need to experiment with the RDP classifier in dada2.

  Citation:

  Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2023): UNITE general FASTA release for Fungi 2. UNITE Community. [DOI:10.15156/BIO/2938068](https://doi.org/10.15156/BIO/2938068)
