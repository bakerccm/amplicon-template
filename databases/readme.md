# Taxonomic Databases

See [scripts/download_databases.sh](../scripts/download_databases.sh) for code to download databases.

General information on training sets for use with dada2: https://benjjneb.github.io/dada2/training.html

## Silva 138.2 training data formatted for DADA2

See [this page](https://zenodo.org/records/14169026) for notes and download links.

  - Use with 16S & assignTaxonomy(): [silva_nr99_v138.2_toGenus_trainset.fa.gz](https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz) (md5:1764e2a36b4500ccb1c7d5261948a414)

  The alternative [silva_nr99_v138.2_toSpecies_trainset.fa.gz](https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toSpecies_trainset.fa.gz) (md5:e0455a1a7de820039d684aad2053e937) is the same sequences but with taxonomy going down to species. But we are assigning species-level taxonomy with addSpecies()/assignSpecies() instead.

  - Use with 16S & assignSpecies(): [silva_v138.2_assignSpecies.fa.gz](https://zenodo.org/records/14169026/files/silva_v138.2_assignSpecies.fa.gz) (md5:62fd939bd1aa9832d2089847f9f0b87d)

  Citations:

  Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucl. Acids Res. 41 (D1): D590-D596. [DOI:10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)

  Yilmaz P, Parfrey LW, Yarza P, Gerken J, Pruesse E, Quast C, Schweer T, Peplies J, Ludwig W, Glöckner FO (2014) The SILVA and "All-species Living Tree Project (LTP)" taxonomic frameworks. Nucl. Acids Res. 42:D643-D648 [DOI:10.1093/nar/gkt1209](https://doi.org/10.1093/nar/gkt1209)

  Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. 2016. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13:581–583. [DOI:10.1038/nmeth.3869](http://doi.org/10.1038/nmeth.3869)
  
  Callahan, BJ (2024). Silva taxonomic training data formatted for DADA2 (Silva version 138.2) [Data set]. Zenodo. [DOI:10.5281/zenodo.14169026](http://doi.org/10.5281/zenodo.14169026)

## UNITE ITS database (Version 10.0 2024-04-04)

See [this page](https://unite.ut.ee/repository.php) for notes and download links (use general FASTA releases with DADA2)

  Fungi, includes global and 97% singletons: [https://doi.org/10.15156/BIO/2959333](https://doi.org/10.15156/BIO/2959333)

  - Use with ITS and assignTaxonomy(): [download](https://s3.hpc.ut.ee/plutof-public/original/db5b0819-2766-4d2f-a2cb-9f100eb447cd.tgz) (md5: a498a463cee353f3faeea00de4447404) This downloads as a tar archive (.tgz) which extracts to two files: sh_general_release_dynamic_s_04.04.2024_dev.fasta (md5: 146699f03719712df7eb2c6bb8358179) and sh_general_release_dynamic_s_04.04.2024.fasta (md5: d43ee757d45db5cb25e7df8d9a08e579). The 'dev' version apparently contains sequence flanking ITS which is removed in the other fasta file. Some classifiers (e.g. the old QIIME2 classifier) are better trained on the full length sequences. May need to experiment with the RDP classifier in dada2.

  Citation:

  Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2024): UNITE general FASTA release for Fungi 2. UNITE Community. [DOI:10.15156/BIO/2959333](https://doi.org/10.15156/BIO/2959333)
