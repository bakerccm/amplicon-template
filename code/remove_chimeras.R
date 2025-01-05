# remove chimeras from sequence table

# load packages

    library("dada2")

# get file paths and parameters from command line arguments

    args = commandArgs(trailingOnly=TRUE)

    filepaths <- list(
        input_rds = args[1], # rds file containing sequence table
        output_rds = args[2]  # rds file containing sequence table without chimeras
    )

# read in sequence table

    sequence.table <- readRDS(file = filepaths$input_rds)

# remove chimeras and save to file

    sequence.table.nochimeras <- removeBimeraDenovo(sequence.table, method = "consensus", multithread = TRUE, verbose = TRUE)

    saveRDS(sequence.table.nochimeras, file = filepaths$output_rds)

# report distribution of sequence lengths

    table(nchar(getSequences(sequence.table.nochimeras)))
