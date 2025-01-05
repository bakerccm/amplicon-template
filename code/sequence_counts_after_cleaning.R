# count sequences post phyloseq

# load packages

    library("phyloseq")

# get file names from command line arguments

    args = commandArgs(trailingOnly=TRUE)

    filenames <- list(
        input = args[1],  # phyloseq object as .rds file
        output = args[2]  # .txt file for output 
    )

    dataset.name <- args[3]
    stage.name <- args[4]

# get input data

    ps <- readRDS(filenames$input)

# count reads and write to output file

    read.sums <- sample_sums(ps)
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
