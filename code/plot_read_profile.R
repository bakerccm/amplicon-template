# plot read quality profiles as PDF

# load packages

    library("dada2")

# get file names from command line arguments

    args = commandArgs(trailingOnly=TRUE)

    filenames <- list(
        input.read1 = args[1],
        input.read2 = args[2],
        output.pdf = args[3]
    )

# plot read quality profiles

    pdf(filenames$output.pdf, width = 8, height = 5)
    plotQualityProfile(c(filenames$input.read1, filenames$input.read2))
    dev.off()
