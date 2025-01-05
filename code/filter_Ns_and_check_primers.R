# filter Ns from data and check primer orientation

# load packages

    library("dada2")
    library("ShortRead")
    library("Biostrings")

# get file names and primers from command line arguments

    args = commandArgs(trailingOnly=TRUE)
    
    filenames <- list()
    primers <- list()
        
    filenames$input.read1 <- args[1]
    filenames$input.read2 <- args[2]
    primers$forward <- args[3]
    primers$reverse <- args[4]
    filenames$output.read1 <- args[5]
    filenames$output.read2 <- args[6]
    filenames$output.summary <- args[7]

# get all orientations of primers
    
    allOrients <- function(primer) {
        # Create all orientations of the input sequence
        require(Biostrings)
        dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
        orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
                     RevComp = reverseComplement(dna))
        return(sapply(orients, toString))  # Convert back to character vector
    }
    
    FWD.orients <- allOrients(primers$forward)
    REV.orients <- allOrients(primers$reverse)

# filter and trim input files
# saves output to files filenames$output.read1 and filenames$output.read2
    
    filterAndTrim(filenames$input.read1, filenames$output.read1, filenames$input.read2, filenames$output.read2, maxN = 0, multithread = TRUE)

# assess results
    primerHits <- function(primer, fn) {
        # Counts number of reads in which the primer is found
        nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
        return(sum(nhits > 0))
    }

    sink(filenames$output.summary, append = FALSE)
    rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = filenames$output.read1),
          FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = filenames$output.read2),
          REV.ForwardReads = sapply(REV.orients, primerHits, fn = filenames$output.read1),
          REV.ReverseReads = sapply(REV.orients, primerHits, fn = filenames$output.read2))
    sink()
