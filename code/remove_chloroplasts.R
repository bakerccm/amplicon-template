# code to filter undesired taxa from 16S dataset

# load packages

    library("here")
    library("phyloseq")
    library("Biostrings")
    library("tidyverse")

# get file names

    args = commandArgs(trailingOnly=TRUE)

    # for debugging
    # args = c(here("out", "16S", "phyloseq", "phyloseq.rds"), here("out", "16S", "phyloseq", "phyloseq_cleaned.rds"), here("out", "16S", "phyloseq", "phyloseq_cleaned.fasta"))

    filenames <- list(
        input_phyloseq = args[1], # .rds file containing 16S phyloseq object
        output_phyloseq = args[2], # .rds file for 16S phyloseq object with chloroplasts, mitochondria etc filtered out
        output_fasta = args[3] # .fasta file for 16S sequences with chloroplasts, mitochondria etc filtered out
    )

# get input data

    phy.16S <- readRDS(filenames$input_phyloseq)

# examine taxonomic assignments

    # extract taxon data as tibbles
    taxa.16S <- tax_table(phy.16S) %>% as("matrix") %>% as_tibble()

    taxa.16S %>% colnames()
    #  [1] "silva_kingdom"          "silva_phylum"           "silva_class"            "silva_order"            "silva_family"           "silva_genus"            "silva_species_single"   "silva_species_multiple" "rdp_kingdom"
    # [10] "rdp_phylum"             "rdp_class"              "rdp_order"              "rdp_family"             "rdp_genus"              "rdp_species_single"     "rdp_species_multiple"   "decipher_domain"        "decipher_phylum"
    # [19] "decipher_class"         "decipher_order"         "decipher_family"        "decipher_genus"         "decipher_species"

    # Eukaryotes
    taxa.16S %>% filter(silva_kingdom == "Eukaryota") # 21 rows; unclear what these are; some might be fungi
    taxa.16S %>% filter(rdp_kingdom == "Eukaryota") # 2 rows of mitochondria; not among those identified as Eukaryota by silva_kingdom, but included among those identified as Mitochondria by silva_kingdom
    taxa.16S %>% filter(decipher_domain == "Eukaryota") # 11 rows; all among those identified by silva_kingdom; unclear what these are; some might be fungi

    # Archaea
    taxa.16S %>% filter(silva_kingdom == "Archaea") # 59 rows
    taxa.16S %>% filter(rdp_kingdom == "Archaea") # 59 rows; same rows as for silva_kingdom
    taxa.16S %>% filter(decipher_domain == "Archaea") # 44 rows, all included within those identified by silva_kingdom and rdp_kingdom

    # mitochondria
    taxa.16S %>% filter(silva_family == "Mitochondria") # 424 rows
    taxa.16S %>% filter(rdp_kingdom == "Mitochondria") # None
    taxa.16S %>% filter(decipher_family == "Mitochondria") # 31 rows; 30 of which included among those identified as Mitochondria by silva_family

    # NAs
    taxa.16S %>% filter(is.na(silva_kingdom)) # None
    taxa.16S %>% filter(is.na(rdp_kingdom)) # 66 rows; 53 of which included among those identified as Mitochondria by silva_kingdom
    taxa.16S %>% filter(is.na(decipher_domain)) # 1,160 rows, many of which seem adequately identified by silva or rdp

    # chloroplasts
    # taken together these filters identify 43 rows as chloroplasts
    taxa.16S %>% filter(silva_order == "Chloroplast" | silva_order == "Cyanobacteria/Chloroplast") # 26 rows
    taxa.16S %>% filter(rdp_phylum == "Chloroplast" | rdp_phylum == "Cyanobacteria/Chloroplast") # 37 rows
    taxa.16S %>% filter(rdp_class == "Chloroplast" | rdp_class == "Cyanobacteria/Chloroplast") # 17 rows
    taxa.16S %>% filter(rdp_order == "Chloroplast" | rdp_order == "Cyanobacteria/Chloroplast") # 19 rows
    taxa.16S %>% filter(decipher_order == "Chloroplast" | decipher_order == "Cyanobacteria/Chloroplast") # 12 rows

# filter out undesired taxa

    ## note that subset_taxa removes NAs ##

    phy.16S.cleaned <- phy.16S %>%
        # Eukaryotes and NAs
            subset_taxa(silva_kingdom != "Eukaryota") %>% # would also remove NAs in silva_kingdom, but there aren't any
            subset_taxa(rdp_kingdom != "Eukaryota") %>%  # also removes NAs in rdp_kingdom
            subset_taxa(decipher_domain != "Eukaryota" | is.na(decipher_domain)) %>% # many of the decipher_domain NAs seem to be assigned just fine with silva or RDP
        # uncomment to remove Archaea
        #    subset_taxa(silva_kingdom != "Archaea" | is.na(silva_kingdom)) %>%
        #    subset_taxa(rdp_kingdom != "Archaea" | is.na(rdp_kingdom)) %>%
        #    subset_taxa(decipher_domain != "Archaea" | is.na(decipher_domain)) %>%
        # mitochondria
            subset_taxa(silva_family != "Mitochondria" | is.na(silva_family)) %>%
            subset_taxa(decipher_family != "Mitochondria" | is.na(decipher_family)) %>%
        # chloroplasts
            subset_taxa(silva_order != "Chloroplast" | is.na(silva_order)) %>%
            subset_taxa(silva_order != "Cyanobacteria/Chloroplast" | is.na(silva_order)) %>%
            subset_taxa(rdp_phylum != "Chloroplast" | is.na(rdp_phylum)) %>%
            subset_taxa(rdp_phylum != "Cyanobacteria/Chloroplast" | is.na(rdp_phylum)) %>%
            subset_taxa(rdp_class != "Chloroplast" | is.na(rdp_class)) %>%
            subset_taxa(rdp_class != "Cyanobacteria/Chloroplast" | is.na(rdp_class)) %>%
            subset_taxa(rdp_order != "Chloroplast" | is.na(rdp_order)) %>%
            subset_taxa(rdp_order != "Cyanobacteria/Chloroplast" | is.na(rdp_order)) %>%
            subset_taxa(decipher_order != "Chloroplast" | is.na(decipher_order)) %>%
            subset_taxa(decipher_order != "Cyanobacteria/Chloroplast" | is.na(decipher_order))

# look at lengths of merged 16S reads (read1 and read2 were previously merged with mergePairs() in script dada.R)

    # in the unfiltered dataset

        refseq(phy.16S) %>% lengths() %>% table() # modal length is 253 bp, consistent with expectations, but there is a tail of shorter and longer sequences

        original_longseqs <- prune_taxa((refseq(phy.16S) %>% lengths()) >= 260, phy.16S)
        tax_table(original_longseqs) # lots of them are mitochondria

        original_shortseqs <- prune_taxa((refseq(phy.16S) %>% lengths()) <= 250, phy.16S)
        tax_table(original_shortseqs) # lots of them are mitochondria

    # in the filtered dataset

        refseq(phy.16S.cleaned) %>% lengths() %>% table() # filtering doesn't totally clean up the length variation in 16S but it is pretty good

        filtered_longseqs <- prune_taxa((refseq(phy.16S.cleaned) %>% lengths()) >= 260, phy.16S.cleaned)
        tax_table(filtered_longseqs) # mostly NAs in here now

        filtered_shortseqs <- prune_taxa((refseq(phy.16S.cleaned) %>% lengths()) <= 250, phy.16S.cleaned)
        tax_table(filtered_shortseqs) # lots of Abditibacteriales in here. This looks like part of candidate phylum FBP, first isolated from ice-free Antarctic soil samples (see https://doi.org/10.1016/j.syapm.2018.01.009)

        # filtered_longseqs %>% refseq() %>% writeXStringSet(filepath = "long_16S_asvs.fasta")

# filter 16S taxa based on fragment lengths

    phy.16S.cleaned <- prune_taxa((refseq(phy.16S.cleaned) %>% lengths()) %in% 250:260, phy.16S.cleaned)

# save cleaned and filtered phyloseq object

    saveRDS(phy.16S.cleaned, file = filenames$output_phyloseq)

# save sequences to FASTA file

    writeXStringSet(refseq(phy.16S.cleaned), filepath = filenames$output_fasta)
