# assign taxonomy using DECIPHER::IdTaxa()

# load packages

    library("dada2")
    library("DECIPHER")

# get file paths from command line arguments

    args = commandArgs(trailingOnly=TRUE)

    filepaths <- list(
        input_file = args[1], # .rds file for input sequence table without chimeras
        decipher_ref = args[2], # .fasta.gz file for DECIPHER::IdTaxa()
        output_file = args[3] # .rds file for output sequence table with taxonomy attached
    )

# read in data

    # sequence table without taxonomy
    sequence.table.nochimeras <- readRDS(filepaths$input_file)

    # DECIPHER reference database
    load(filepaths$decipher_ref) # loads as trainingSet

# get sequences from ASV table and convert to DNAStringSet

    dna <- DNAStringSet(getSequences(sequence.table.nochimeras))

# assign taxonomy with DECIPHER

    ids <- IdTaxa(dna, trainingSet, strand = "top", processors = NULL, verbose = FALSE) # use all processors

# convert the DECIPHER output (class "Taxa") to a matrix analogous to the output from dada2::assignTaxonomy()

    ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
    taxid <- t(sapply(ids, function(x) {
            m <- match(ranks, x$rank)
            taxa <- x$taxon[m]
            taxa[startsWith(taxa, "unclassified_")] <- NA
            taxa
    }))
    colnames(taxid) <- ranks
    rownames(taxid) <- getSequences(sequence.table.nochimeras)

    # DADA2: "The taxid matrix from IdTaxa is a drop-in replacement for the taxa matrix from assignTaxonomy, simply set taxa <- taxid to carry on using the IdTaxa assignments."

# save taxa to file

    saveRDS(taxid, file = filepaths$output_file)
