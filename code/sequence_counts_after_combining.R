# count sequences post combining (should be the same as post-chimeras)

# load packages

    library("phyloseq")

# get file names from command line arguments

    args = commandArgs(trailingOnly=TRUE)

    filenames <- list(
        input = args[1],  # list of phyloseq objects as .rds file
        output = args[2]  # csv file for output 
    )

    stage.name <- args[3]

# get input data

    ps.list <- readRDS(filenames$input) # should contain list of phyloseq objects

# count reads and write to output file

read.sums.list <- list()

for (dataset in names(ps.list)) {

    read.sums <- sample_sums(ps.list[[dataset]])
    read.sums.list[[dataset]] <- data.frame(
        dataset = dataset,
        stage = stage.name,
        sample = names(read.sums),
        reads = read.sums
    )
    rownames(read.sums.list[[dataset]]) <- NULL

}

read.sums.df <- do.call("rbind", read.sums.list)

write.table(
    read.sums.df,
    file = filenames$output,
    sep = ",", row.names = FALSE, quote= FALSE
)

