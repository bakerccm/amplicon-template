# sample inference using dada function of dada2

# load packages

    library("dada2")

# get file paths and parameters from command line arguments

    args = commandArgs(trailingOnly=TRUE)

    filepaths <- list(
        input_fastqs = list(
            read1 = args[1], # fastq.gz file for read 1
            read2 = args[2]  # fastq.gz file for read 2
        ),
        input_error_profiles = list(
            read1 = args[3], # rds file containing estimated error profile appropriate for input_fastqs$read1 (i.e. same sequencing run and read direction)
            read2 = args[4]  # rds file containing estimated error profile appropriate for input_fastqs$read2 (i.e. same sequencing run and read direction)
        ),
        output_rds = args[5]  # save mergePairs() output to rds file
    )

# read in previously-computed error profiles

    error_profile <- sapply(filepaths$input_error_profiles, simplify=FALSE, readRDS)

# sample inference

    # 5 seconds
    derep <- sapply(filepaths$input_fastqs, simplify=FALSE, derepFastq)

    # 30 seconds
    dd <- sapply(names(derep), simplify=FALSE, function (X) dada(derep = derep[[X]], err = error_profile[[X]], multithread = TRUE))

# merge paired end reads

    # 2 seconds
    merged.pairs <- mergePairs(dadaF = dd$read1, derepF = derep$read1, dadaR = dd$read2, derepR = derep$read2)

# save merged pairs

    saveRDS(merged.pairs, filepaths$output_rds)
