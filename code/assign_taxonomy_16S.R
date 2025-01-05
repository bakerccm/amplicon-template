# assign taxonomy using assignTaxonomy() and addSpecies() [which wraps assignSpecies()]
# --> use for 16S dataset

# load packages

    library("dada2")

# get file paths from command line arguments

    args = commandArgs(trailingOnly=TRUE)

    filepaths <- list(
        input_file = args[1], # .rds file for input sequence table without chimeras
        assignTaxonomy_refFasta = args[2], # .fasta.gz file for assignTaxonomy()
        assignSpecies_refFasta = args[3], # .fasta.gz file for assignSpecies()
        output_file_single = args[4], # .rds file for output sequence table with taxonomy attached
        output_file_multiple = args[5] # .rds file for output sequence table with taxonomy attached (multiple matches allowed)
    )

# read in data

    sequence.table.nochimeras <- readRDS(filepaths$input_file)

# assign taxonomy to genus level
#   minBoot = 50 by default (recommended for sequences of 250bp and less; see Wang et al 2007)
#   tryRC = FALSE by defaultl; consider using TRUE if many sequences fail to assign

    taxa <- assignTaxonomy(seqs = sequence.table.nochimeras, refFasta = filepaths$assignTaxonomy_refFasta, multithread = TRUE)

# add species assignments where genus is consistent
# tryRC = FALSE by default; consider using TRUE if many sequences fail to assign

    # only allow unambiguous identifications
    taxa.species.single <- addSpecies(taxtab = taxa, refFasta = filepaths$assignSpecies_refFasta, allowMultiple = FALSE)

    # allow multiple species matches
    taxa.species.multiple <- addSpecies(taxtab = taxa, refFasta = filepaths$assignSpecies_refFasta, allowMultiple = TRUE)

# save taxa to file

    saveRDS(taxa.species.single, file = filepaths$output_file_single)

    saveRDS(taxa.species.multiple, file = filepaths$output_file_multiple)
