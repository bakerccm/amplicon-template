# code to combine data incl 16S SEPP tree into a single phyloseq object

# load packages

    library("tidyverse")

    # note phyloseq install instructions at https://joey711.github.io/phyloseq/install.html use BiocLite rather than BiocManager as used here; the latter is preferred with R > 3.5
    # library("BiocManager")
    # BiocManager::install("phyloseq")
    library("phyloseq")

# get file names

    args = commandArgs(trailingOnly=TRUE)

    # for debugging
    # args <- c(
    #     file.path("out", "16S", "phyloseq", "phyloseq_cleaned.rds"),
    #     file.path("out", "16S", "phyloseq", "phyloseq_cleaned_placement.tog.tre"),
    #     file.path("out", "ITS", "out", "ITS", "phyloseq", "phyloseq.rds"),
    #     file.path("out", "combined", "amplicon.rds")
    # )

    filenames <- list(
        "16S" = args[1], # .rds file for 16S phyloseq object with chloroplasts, mitochondria etc filtered out
        "16S_sepp_tree" = args[2], # 16S placement tree from SEPP
        "ITS" = args[3], # .rds file containing ITS phyloseq object
        "output" = args[4] # .rds file containing original and normalized 16S and ITS phyloseq objects
    )

# get data, store as list of phyloseq objects

    ps.list <- lapply(filenames[c("16S", "ITS")], readRDS)

# attach SEPP tree for 16S data

    sepp_tree_16S <- read_tree_greengenes(filenames$`16S_sepp_tree`)

    ps.list$`16S` <- merge_phyloseq(ps.list$`16S`, phy_tree(sepp_tree_16S))

# save ps.list i.e. list of phyloseq objects

    saveRDS(ps.list, file = filenames$output)
