# code to filter undesired taxa from 16S dataset

# load packages

    library("phyloseq")
    library("Biostrings")

# get file names

    args = commandArgs(trailingOnly=TRUE)

    filenames <- list(
        input_phyloseq = args[1], # .rds file containing 16S phyloseq object
        output_phyloseq = args[2], # .rds file for 16S phyloseq object with chloroplasts, mitochondria etc filtered out
        output_fasta = args[3] # .fasta file for 16S sequences with chloroplasts, mitochondria etc filtered out
    )

# get input data

    phy.16S <- readRDS(filenames$input_phyloseq)

# filter out undesired taxa

    ## note that subset_taxa removes NAs ##

    # Eukaryotes and NAs
    phy.16S.cleaned <- subset_taxa(phy.16S, silva_kingdom != "Eukaryota")

    # uncomment to remove Archaea
    #   phy.16S.cleaned <- subset_taxa(phy.16S.cleaned, silva_kingdom != "Archaea" | is.na(silva_kingdom))
   
    # mitochondria
    phy.16S.cleaned <- subset_taxa(phy.16S.cleaned, silva_family != "Mitochondria" | is.na(silva_family))
    
    # chloroplasts
    phy.16S.cleaned <- subset_taxa(phy.16S.cleaned, silva_order != "Chloroplast" | is.na(silva_order))
    phy.16S.cleaned <- subset_taxa(phy.16S.cleaned, silva_order != "Cyanobacteria/Chloroplast" | is.na(silva_order))

# filter 16S taxa based on fragment lengths

    phy.16S.cleaned <- prune_taxa(lengths(refseq(phy.16S.cleaned)) %in% 250:260, phy.16S.cleaned)

# save cleaned and filtered phyloseq object

    saveRDS(phy.16S.cleaned, file = filenames$output_phyloseq)

# save sequences to FASTA file

    writeXStringSet(refseq(phy.16S.cleaned), filepath = filenames$output_fasta)
