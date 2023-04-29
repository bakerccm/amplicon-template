# count sequences post normalization

# load packages

    library("phyloseq")

# get file names from command line arguments

    args = commandArgs(trailingOnly=TRUE)

    filenames <- list(
        input = args[1],  # list of phyloseq objects as .rds file
        output.dir = args[2]  # directory for output files 
    )

# get input data

    load(filenames$input) # should contain ps.list and ps.list.tss

# count reads and write to output file

for dataset in names(ps.list) {

    read.sums <- sample_sums(ps.list[[dataset]])
    read.sums.df <- data.frame(sample = names(read.sums), reads = read.sums)
    rownames(read.sums.df) <- NULL
    current.filename <- paste0(dataset, ".txt")
    write.table(read.sums.df, file = paste(filenames$output.dir, current.filename, sep = "/"), row.names = FALSE, quote= FALSE)

}
