# filter 16S datasets
# uses maxN, truncQ, maxEE, truncLen arguments

# load packages

    library("dada2")

# get file names from command line arguments

    args = commandArgs(trailingOnly=TRUE)

    filenames <- list(
        input.read1 = args[1],
        input.read2 = args[2],
        output.read1 = args[3],
        output.read2 = args[4]
    )

    params <- list(
        maxN = as.numeric(args[5]),
        truncQ = as.numeric(args[6]),
        maxEE = c(as.numeric(args[7]), as.numeric(args[8])),
        truncLen = c(as.numeric(args[9]), as.numeric(args[10]))
    )

# filter reads

    filterAndTrim(filenames$input.read1, filenames$output.read1, filenames$input.read2, filenames$output.read2,
        maxN = params$maxN, truncQ = params$truncQ, maxEE = params$maxEE, truncLen = params$truncLen,
        rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
