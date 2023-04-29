# count sequences in sequence table (post dada, post chimeras etc)

# load packages

    library("dada2")

# get file names from command line arguments

    args = commandArgs(trailingOnly=TRUE)

    filenames <- list(
        input = args[1],  # sequence table as .rds file
        output = args[2]  # .txt file for output 
    )

# get input data

    sequence.table <- readRDS(filenames$input)

# count reads and write to output file

    read.sums <- apply(sequence.table, MAR = 1 , sum)
    read.sums.df <- data.frame(sample = names(read.sums), reads = read.sums)
    rownames(read.sums.df) <- NULL
    write.table(read.sums.df, file = filenames$output, row.names = FALSE, quote = FALSE)
