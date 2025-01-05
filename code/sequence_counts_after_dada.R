# count sequences in sequence table (post dada, post chimeras etc)

# load packages

    library("dada2")

# get file names from command line arguments

    args = commandArgs(trailingOnly=TRUE)

    filenames <- list(
        input = args[1],  # sequence table as .rds file
        output = args[2]  # .csv file for output 
    )

    dataset.name <- args[3]
    stage.name <- args[4]

# get input data

    sequence.table <- readRDS(filenames$input)

# count reads and write to output file

    read.sums <- apply(sequence.table, MAR = 1 , sum)
    read.sums.df <- data.frame(
        dataset = dataset.name,
        stage = stage.name,
        sample = names(read.sums),
        reads = read.sums
    )
    rownames(read.sums.df) <- NULL
    write.table(
        read.sums.df,
        file = filenames$output,
        sep = ",", row.names = FALSE, quote = FALSE
    )
