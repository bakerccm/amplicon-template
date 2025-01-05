# learn error profile for files in directory

# input directory should contain either read1 or read2 (not both) fastq.gz files from a single run (i.e. files with the same error profile)

# note that learnErrors() starts with the first file in the directory, and continues to read in files in order until 'nbases' bases are
# obtained (set in config file; default is 1e8); thus file order matters, though hopefully estimated error profile is similar regardless; files can
# be read in random order using the randomize option, but this still only samples from a subset of files (i.e. does not subsample from all files)

# load packages

    library("dada2")

# get file paths and parameters from command line arguments

    args = commandArgs(trailingOnly=TRUE)

    filepaths <- list(
        input_file_dir = args[1],
        output_rds = args[2], # save error profile as rds file
        output_pdf = args[3]  # save PDF of estimated error rates
    )

    params <- list(
        nbases = as.numeric(args[4])
    )

# estimate error profile

    estimated_error_profile <- learnErrors(filepaths$input_file_dir, nbases = params$nbases, multithread = TRUE)

# save error profile to file

    saveRDS(estimated_error_profile, file = filepaths$output_rds)

# plot errors

    pdf(filepaths$output_pdf, height = 6, width = 6)
    plotErrors(estimated_error_profile, nominalQ = TRUE)
    dev.off()
