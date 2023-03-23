# Amplicon data project

## Data description

### Overview

This amplicon sequencing project comprises DNA sequence reads generated from <experimental details to be supplied here>.

DNA extracted from each sample was sent to the Environmental Sample Preparation and Sequencing Facility at Argonne National Laboratory (ANL; Lemont IL 60439, United States) for PCR amplification, multiplexed library preparation and Illumina MiSeq sequencing. ANL used primers targeting the 16S rRNA gene and the Internal Transcribed Spacer (ITS) genetic region in order to profile bacterial and fungal communities, respectively. Key details and citations are provided below.

The pipeline in this repo is designed to process the six files (forward, reverse and index reads for the 16S and ITS datasets) provided by ANL to CRREL prior to any demultiplexing or analysis.

### 16S dataset

The 16S dataset was generated following the protocol outlined in [Caporaso *et al*. 2011](#citations), with minor modifications to some quantities and to the primers used:

**515FY:** 5'-GTGYCAGCMGCCGCGGTAA-3' ([Parada *et al*. 2016](#citations); sometimes referred to as 515F)

**806RB:** 5'-GGACTACNVGGGTWTCTAAT-3' ([Apprill *et al*. 2015](#citations); sometimes referred to as 806R)

These 515FY ([Parada *et al*. 2016](#citations)) and 806RB ([Apprill *et al*. 2015](#citations)) primers are modified versions of the original F515 (5'-GTGCCAGCMGCCGCGGTAA-3') and R806 (5'-GGACTACHVGGGTWTCTAAT-3') primers from [Caporaso *et al*. 2011](#citations). Forward-barcoded constructs also incorporate changes described by [Walters *et al*. 2016](#citations) based upon those in [Caporaso *et al*. 2012](#citations). According to Caporaso *et al*., these primers "target the V4 region of the 16S rRNA and, for reference, amplify the region 533–786 in *E.coli* strain 83972 (greengenes accession no. prokM- SA_id:470367)". Apprill *et al*. modified Caporaso *et al*.'s original primer with the goal of enhancing capture of SAR11 targets. Parada *et al*. pair 515FY with a different primer 926R (5'-CCGYCAATTYMTTTRAGTTT-3'); however, according to Argonne (24 May 2021) in their experience many researchers have not observed the improved resolution reported by Parada *et al*., and Argonne still most frequently uses 515FY/806RB.

### ITS dataset

The ITS dataset was generated following the protocol outlined in Smith and Peay 2014, with minor modification to the some quantities. This protocol used ITS1F (Gardes and Bruns 1993) and ITS2 primers ([White *et al* 1990](#citations)):

**ITS1F:** 5'-CTTGGTCATTTAGAGGAAGTAA-3' ([Gardes and Bruns 1993](#citations))

**ITS2:** 5'-GCTGCGTTCTTCATCGATGC-3' ([White *et al* 1990](#citations))

These primers target the ITS1 region from the fungal nuclear ribosomal RNA gene. [Smith and Peay 2014](#citations) used them as their 'PCR primers', leaving the sequences unmodified, but adding Illumina adaptors etc on the 5' ends. Note, however, that Smith and Peay used additional 'sequencing primers' that overlap and extend beyond the 3' ends of the PCR primers:

**Read1 seq primer:** 5'-TTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC-3' ([Smith and Peay 2014](#citations))

**Read2 seq primer:** 5'-CGTTCTTCATCGATGCVAGARCCAAGAGATC-3' ([Smith and Peay 2014](#citations))

### Sequencing

Illumina Miseq 2 x 251bp paired end + 12bp index

### Citations

Apprill, Amy, Sean McNally, Rachel Parsons and Laura Weber (2015). Minor revision to V4 region SSU rRNA 806R gene primer greatly increases detection of SAR11 bacterioplankton. *Aquatic Microbial Ecology* 75:129-137. doi:10.3354/ame01753

Caporaso, J Gregory, Christian L Lauber, William A Walters, Donna Berg-Lyons, Catherine A Lozupone, Peter J Turnbaugh, Noah Fierer and Rob Knight (2011). Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample. *PNAS* 108(Supplement 1):4516-4522. doi:10.1073/pnas.1000080107

Caporaso, J Gregory, Christian L Lauber, William A Walters, Donna Berg-Lyons, James Huntley, Noah Fierer, Sarah M Owens, Jason Betley, Louise Fraser, Markus Bauer, Niall Gormley, Jack A Gilbert, Geoff Smith and Rob Knight (2012). Ultra-high-throughput microbial community analysis on the Illumina HiSeq and MiSeq platforms. *The ISME journal* 6(8):1621-1624. doi:10.1038/ismej.2012.8

Gardes, M and T Bruns (1993). ITS primers with enhanced specificity for basidiomycetes - application to the identification of mycorrhizae and rusts. *Molecular Ecology* 2:113–118. doi:10.1111/j.1365-294x.1993.tb00005.x

Parada, Alma E, David M Needham and Jed A Fuhrman (2016). Every base matters: assessing small subunit rRNA primers for marine microbiomes with mock communities, time series and global field samples. *Environmental Microbiology* 18:1403-1414. doi:10.1111/1462-2920.13023

Smith, Dylan P and Kabir G Peay (2014). Sequence Depth, Not PCR Replication, Improves Ecological Inference from Next Generation DNA Sequencing. *PLoS ONE* 9(2):e90234. doi:10.1371/journal.pone.0090234

Walters, William, Embriette R Hyde, Donna Berg-Lyons, Gail Ackermann, Greg Humphrey, Alma Parada, Jack A Gilbert, Janet K Jansson, J Gregory Caporaso, Jed A Fuhrman, Amy Apprill and Rob Knight (2016). Improved Bacterial 16S rRNA Gene (V4 and V4-5) and Fungal Internal Transcribed Spacer Marker Gene Primers for Microbial Community Surveys. *mSystems* 1(1):e00009-15. doi:10.1128/mSystems.00009-15

White, TJ, TD Bruns, S Lee and J Taylor (1990). Amplification and direct sequencing of fungal ribosomal RNA genes for phylogenetics. In: Innis MA, DH Gefland, JJ Sninsky and TJ White (eds). *PCR protocols: a guide to method and applications*. San Diego, Academic Press. pp. 315-322. doi:10.1016/b978-0-12-372180-8.50042-1

---
Created by [Chris Baker](https://github.com/bakerccm)
