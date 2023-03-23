# package sequence data up with sample and taxonomy metadata as a phyloseq object ready for analysis
# -- use for 16S dataset

# load packages

    library("phyloseq")
    library("Biostrings")

# get file names from command line arguments

    args = commandArgs(trailingOnly=TRUE)

    filenames <- list(
        input_file_sequence_table = args[1], # .rds file for sequence table without chimeras
        input_file_silva_single = args[2], # .rds file for silva taxonomy assignments (unambiguous species assignments only)
        input_file_silva_multiple = args[3], # .rds file for silva taxonomy assignments (multiple species assignments allowed)
        input_file_rdp_single = args[4], # .rds file for rdp taxonomy assignments (unambiguous species assignments only)
        input_file_rdp_multiple = args[5], # .rds file for rdp taxonomy assignments (multiple species assignments allowed)
        input_file_decipher = args[6], # .rds file for DECIPHER taxonomy assignments
        input_file_metadata = args[7], # sample metadata file from config/
        output_phyloseq = args[8], # .rds file name to save phyloseq object output
        output_fasta = args[9] # .fasta file name to save sequences to
    )

# get sequence table

    sequence.table <- readRDS(filenames$input_file_sequence_table)

    # replace dashes with periods in row names to avoid headaches in phyloseq later
        rownames(sequence.table) <- gsub("-", ".", rownames(sequence.table))

# get taxonomy results

    taxa <- list(
        silva_single = readRDS(filenames$input_file_silva_single),
        silva_multiple = readRDS(filenames$input_file_silva_multiple),
        rdp_single = readRDS(filenames$input_file_rdp_single),
        rdp_multiple = readRDS(filenames$input_file_rdp_multiple),
        decipher = readRDS(filenames$input_file_decipher)
    )

# sanity checks

    # single- and multiple-match results should be the same over the Kingdom:Genus fields
        try(if(!identical(taxa$silva_single[,1:6], taxa$silva_multiple[,1:6])) stop("silva_single and silva_multiple do not match on Kingdom:Genus"))
        try(if(!identical(taxa$rdp_single[,1:6], taxa$rdp_multiple[,1:6])) stop("rdp_single and rdp_multiple do not match on Kingdom:Genus"))

    # columns of sequence.table should match rows of taxonomic assignment tables (though if they don't, it could probably be fixed with a match or a join)
        try(if(!identical(colnames(sequence.table), rownames(taxa$silva_single))) stop("sequence.table column names do not match silva_single row names"))
        try(if(!identical(colnames(sequence.table), rownames(taxa$silva_multiple))) stop("sequence.table column names do not match silva_multiple row names"))
        try(if(!identical(colnames(sequence.table), rownames(taxa$rdp_single))) stop("sequence.table column names do not match rdp_single row names"))
        try(if(!identical(colnames(sequence.table), rownames(taxa$rdp_multiple))) stop("sequence.table column names do not match rdp_multiple row names"))
        try(if(!identical(colnames(sequence.table), rownames(taxa$decipher))) stop("sequence.table column names do not match decipher row names"))

# get sample metadata and parse

    sample.metadata <- read.table(filenames$input_file_metadata, header = TRUE, stringsAsFactors = FALSE)

    # replace dashes with periods to avoid headaches in phyloseq later
        sample.metadata$SampleID <- gsub("-", ".", sample.metadata$SampleID)

    # rownames need to match rownames(sequence.table) when creating phyloseq object later
        rownames(sample.metadata) <- gsub("-", ".", sample.metadata$SampleID)

    # filter to 16S only
        sample.metadata <- sample.metadata[sample.metadata$Dataset == "16S",]

    # remove unnecessary columns
        remove <- c("BarcodeSequence", "LinkerPrimerSequence")
        sample.metadata <- sample.metadata[!(names(sample.metadata) %in% remove)]

    # match rows to sequence table
        sample.metadata <- sample.metadata[base::match(rownames(sequence.table), rownames(sample.metadata)),]

# prefix/suffix taxonomy colnames with taxonomy assignment method and then cbind the different tables together

    colnames(taxa$silva_single) <- c(paste0("silva_", colnames(taxa$silva_single)[1:6]), paste0("silva_", colnames(taxa$silva_single)[7], "_single"))
    colnames(taxa$silva_multiple) <- c(paste0("silva_", colnames(taxa$silva_multiple)[1:6]), paste0("silva_", colnames(taxa$silva_multiple)[7], "_multiple"))

    colnames(taxa$rdp_single) <- c(paste0("rdp_", colnames(taxa$rdp_single)[1:6]), paste0("rdp_", colnames(taxa$rdp_single)[7], "_single"))
    colnames(taxa$rdp_multiple) <- c(paste0("rdp_", colnames(taxa$rdp_multiple)[1:6]), paste0("rdp_", colnames(taxa$rdp_multiple)[7], "_multiple"))

    colnames(taxa$decipher) <- paste0("decipher_", colnames(taxa$decipher))

    taxonomy.metadata <- cbind(taxa$silva_single, taxa$silva_multiple[,7,drop=FALSE], taxa$rdp_single, taxa$rdp_multiple[,7,drop=FALSE], taxa$decipher)

    colnames(taxonomy.metadata) <- tolower(colnames(taxonomy.metadata))

# construct phyloseq object from sequence table, sample metadata and taxonomy metadata

    ps <- phyloseq(otu_table(sequence.table, taxa_are_rows=FALSE),
        sample_data(sample.metadata),
        tax_table(taxonomy.metadata))

# rename ASVs and store DNA sequences in the refseq slot of the phyloseq object [can be accessed later with refseq(ps)]

    dna <- Biostrings::DNAStringSet(taxa_names(ps))
    names(dna) <- taxa_names(ps)
    ps <- merge_phyloseq(ps, dna)
    taxa_names(ps) <- paste0("ASV_16S_", sprintf("%05d",1:ntaxa(ps)))

# save phyloseq object

    saveRDS(ps, file = filenames$output_phyloseq)

# save sequences to FASTA file

    writeXStringSet(refseq(ps), filepath = filenames$output_fasta)
