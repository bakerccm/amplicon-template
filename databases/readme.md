# Taxonomic Databases

General information on training sets for use with dada2: https://benjjneb.github.io/dada2/training.html

## Silva 138.1 prokaryotic SSU taxonomic training data formatted for DADA2

See this page for notes and download links: https://zenodo.org/record/4587955#.YMEPdC1r3UI

  - Use with 16S & assignTaxonomy(): [silva_nr99_v138.1_train_set.fa.gz](https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz) (md5:6b41db7139834c71171f8ce5b5918fc6)

  The alternative [silva_nr99_v138.1_wSpecies_train_set.fa.gz](https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz) (md5:ba13ab369161e0cb85df7e0ee3a4182e) is the same sequences but with taxonomy going down to species. But we are assigning species-level taxonomy with addSpecies()/assignSpecies() instead.

  - Use with 16S & assignSpecies(): [silva_species_assignment_v138.1.fa.gz](https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz) (md5:f21c2d97c79ff07c17949a9622371a4c)

  - Use with 16S & DECIPHER (modified SILVA SSU r138 dataset): [SILVA_SSU_r138_2019.RData](http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData) (md5:cb983b6a5e8cdb46f8c88b5afae21f66)

  Citations:

  Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucl. Acids Res. 41 (D1): D590-D596.

  Yilmaz P, Parfrey LW, Yarza P, Gerken J, Pruesse E, Quast C, Schweer T, Peplies J, Ludwig W, Glöckner FO (2014) The SILVA and "All-species Living Tree Project (LTP)" taxonomic frameworks. Nucl. Acids Res. 42:D643-D648

  Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. 2016. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13:581–583. doi:10.1038/nmeth.3869
  
  McLaren, Michael R., & Callahan, Benjamin J. (2021). Silva 138.1 prokaryotic SSU taxonomic training data formatted for DADA2 [Data set]. Zenodo. http://doi.org/10.5281/zenodo.4587955

## RDP taxonomic training data formatted for DADA2 (RDP trainset 18/release 11.5)

See this page for notes and download links: https://zenodo.org/record/4310151#.YMEX9C1r3UI

  - Use with 16S & assignTaxonomy(): [rdp_train_set_18.fa.gz](https://zenodo.org/record/4310151/files/rdp_train_set_18.fa.gz) (md5:693c1ade0caf223e1fb7c6cf4da6f272)

  - Use with 16S & assignSpecies(): [rdp_species_assignment_18.fa.gz](https://zenodo.org/record/4310151/files/rdp_species_assignment_18.fa.gz) (md5:1a465122bdf670a74eb8873de042bad6)
  
  Citation:
  
  Callahan, Benjamin. (2020). RDP taxonomic training data formatted for DADA2 (RDP trainset 18/release 11.5) [Data set]. Zenodo. http://doi.org/10.5281/zenodo.4310151
  
## UNITE ITS database (Version 8.3 2021-05-10)

See this page for notes and download links: https://unite.ut.ee/repository.php (use general FASTA releases with DADA2)

  Fungi, includes global and 97% singletons. https://doi.org/10.15156/BIO/1280089

  Downloads available from https://plutof.ut.ee/#/doi/10.15156/BIO/1280089

  - Use with ITS and assignTaxonomy(): [download](https://files.plutof.ut.ee/public/orig/6A/F9/6AF94919CCB48307734D6256CACA50AE1ECBC0839F644D4B661E3673525E41A4.tgz) (md5: 51ffc4ab953c2036f72ac5b6cd06f168) This downloads as a tar archive (.tgz) which extracts to a directory `sh_general_release_s_10.05.2021` with two files in it: `sh_general_release_dynamic_s_10.05.2021_dev.fasta` (md5:7f6cdd2dc1b64f115df05abde221b377) and `sh_general_release_dynamic_s_10.05.2021.fasta` (md5:1a3125790488fd74aa5b6ef3d6e9ae79). The 'dev' version apparently contains sequence flanking ITS which is removed in the other fasta file. Some classifiers (e.g. the old QIIME2 classifier) are better trained on the full length sequences. May need to experiment with the RDP classifier in dada2.

  Citation:

  Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2021): UNITE general FASTA release for Fungi 2. Version 10.05.2021. UNITE Community. https://doi.org/10.15156/BIO/1280089
