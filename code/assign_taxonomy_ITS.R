# assign taxonomy using assignTaxonomy() only
# --> use for ITS dataset

# load packages

    library("dada2")

# get file paths from command line arguments

    args = commandArgs(trailingOnly=TRUE)

    filepaths <- list(
        input_file = args[1], # .rds file for input sequence table without chimeras
        assignTaxonomy_refFasta = args[2], # .fasta.gz file for assignTaxonomy()
        output_file = args[3] # .rds file for taxonomy assignments
    )

# read in data

    sequence.table.nochimeras <- readRDS(filepaths$input_file)

# assign taxonomy to genus level
#   minBoot = 50 by default (recommended for sequences of 250bp and less; see Wang et al 2007)
#   tryRC = FALSE by defaultl; consider using TRUE if many sequences fail to assign

    taxa <- assignTaxonomy(seqs = sequence.table.nochimeras, refFasta = filepaths$assignTaxonomy_refFasta, multithread = TRUE)

# save taxonomy assignments to file (note this does not include read counts, so this needs to be interpreted along with the input sequence table)

    saveRDS(taxa, file = filepaths$output_file)
