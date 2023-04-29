# count sequences post phyloseq

# load packages

    library("phyloseq")

# get file names from command line arguments

    args = commandArgs(trailingOnly=TRUE)

    filenames <- list(
        input = args[1],  # phyloseq object as .rds file
        output = args[2]  # .txt file for output 
    )

# get input data

    ps <- readRDS(filenames$input)

# count reads and write to output file

    read.sums <- sample_sums(ps)
    read.sums.df <- data.frame(sample = names(read.sums), reads = read.sums)
    rownames(read.sums.df) <- NULL
    write.table(read.sums.df, file = filenames$output, row.names = FALSE, quote = FALSE)
