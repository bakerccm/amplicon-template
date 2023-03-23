# package sequence data up with sample and taxonomy metadata as a phyloseq object ready for analysis
# -- use for 16S dataset

# load packages

    library("phyloseq")
    library("Biostrings")

# get file names from command line arguments

    args = commandArgs(trailingOnly=TRUE)

    filenames <- list(
        input_file_sequence_table = args[1], # .rds file for sequence table without chimeras
        input_file_unite = args[2], # .rds file for unite taxonomy assignments
        input_file_metadata = args[3], # sample metadata file from config/
        output_file = args[4] # .rds file name to save phyloseq object output
    )

# get sequence table

    sequence.table <- readRDS(filenames$input_file_sequence_table)

    # replace dashes with periods in row names to avoid headaches in phyloseq later
        rownames(sequence.table) <- gsub("-", ".", rownames(sequence.table))

# get taxonomy results

    taxa <- list(
        unite = readRDS(filenames$input_file_unite)
    )

# sanity checks

    # columns of sequence.table should match rows of taxonomic assignment table (though if they don't, it could probably be fixed with a match or a join)
        try(if(!identical(colnames(sequence.table), rownames(taxa$unite))) stop("sequence.table column names do not match decipher row names"))

# get sample metadata and parse

    sample.metadata <- read.table(filenames$input_file_metadata, header = TRUE, stringsAsFactors = FALSE)

    # replace dashes with periods to avoid headaches in phyloseq later
        sample.metadata$SampleID <- gsub("-", ".", sample.metadata$SampleID)

    # rownames need to match rownames(sequence.table) when creating phyloseq object later
        rownames(sample.metadata) <- gsub("-", ".", sample.metadata$SampleID)

    # filter to ITS only
        sample.metadata <- sample.metadata[sample.metadata$Dataset == "ITS",]

    # remove unnecessary columns
        remove <- c("BarcodeSequence", "LinkerPrimerSequence")
        sample.metadata <- sample.metadata[!(names(sample.metadata) %in% remove)]

    # match rows to sequence table
        sample.metadata <- sample.metadata[base::match(rownames(sequence.table), rownames(sample.metadata)),]

# prefix/suffix taxonomy colnames with taxonomy assignment method and then cbind the different tables together

    colnames(taxa$unite) <- paste0("unite_", colnames(taxa$unite))

    taxonomy.metadata <- taxa$unite

    colnames(taxonomy.metadata) <- tolower(colnames(taxonomy.metadata))

# construct phyloseq object from sequence table, sample metadata and taxonomy metadata

    ps <- phyloseq(otu_table(sequence.table, taxa_are_rows=FALSE),
        sample_data(sample.metadata),
        tax_table(taxonomy.metadata))

# rename ASVs and store DNA sequences in the refseq slot of the phyloseq object [can be accessed later with refseq(ps)]

    dna <- Biostrings::DNAStringSet(taxa_names(ps))
    names(dna) <- taxa_names(ps)
    ps <- merge_phyloseq(ps, dna)
    taxa_names(ps) <- paste0("ASV_ITS_", sprintf("%05d",1:ntaxa(ps)))

# save phyloseq object for later use!

    saveRDS(ps, file = filenames$output_file)
