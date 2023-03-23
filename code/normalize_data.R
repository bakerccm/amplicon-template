# code to normalize reads by total sum scaling, add SEPP tree, and
# package up data from 16S and ITS datasets into a single rds file

# load packages

    library("here")
    library("tidyverse")

    # note phyloseq install instructions at https://joey711.github.io/phyloseq/install.html use BiocLite rather than BiocManager as used here; the latter is preferred with R > 3.5
    # library("BiocManager")
    # BiocManager::install("phyloseq")
    library("phyloseq")

# get file names

    args = commandArgs(trailingOnly=TRUE)

    # for debugging
    # args <- c(
    #     here("out", "16S", "phyloseq", "phyloseq_cleaned.rds"),
    #     here("out", "16S", "phyloseq", "phyloseq_cleaned_placement.tog.tre"),
    #     here("out", "ITS", "phyloseq", "phyloseq.rds"),
    #     here("out", "combined", "amplicon_normalized.rdata")
    # )

    filenames <- list(
        "16S" = args[1], # .rds file for 16S phyloseq object with chloroplasts, mitochondria etc filtered out
        "16S_sepp_tree" = args[2], # 16S placement tree from SEPP
        "ITS" = args[3], # .rds file containing ITS phyloseq object
        "output" = args[4] # .rdata file containing original and normalized 16S and ITS phyloseq objects
    )

# get data, store as list of phyloseq objects

    ps.list <- lapply(filenames[c("16S", "ITS")], readRDS)

# attach SEPP tree for 16S data

    sepp_tree_16S <- read_tree_greengenes(filenames$`16S_sepp_tree`)

    ps.list$`16S` <- merge_phyloseq(ps.list$`16S`, phy_tree(sepp_tree_16S))

# Normalize sequence tables to relative read counts (i.e. total sum scaling = TSS),
# and aggregate to genus, family, etc
# (see, e.g., code at https://github.com/rthanis/animal-cdiff/blob/master/analysis/alpha-beta.Rmd)

# N.B. tax_glom removes ASVs without the relevant taxonomic assignment
# so this code calculates total sum scaling once (at ASV level, so without tax_glom)
# and then uses those relative abundances for the aggregated versions, so the
# latter will likely have relative abundances that do not sum to 1 (per sample).

tss_ASVs <- list(
    `16S_ASVs` = ps.list$`16S` %>% transform_sample_counts(function (x) x / sum(x)),
    `ITS_ASVs` = ps.list$`ITS` %>% transform_sample_counts(function (x) x / sum(x))
)

# do each database separately, since the taxon ranks to the left of the column used for tax_glom should be only the hierarchy for the database being glommed

    tax_table_16S <- tax_table(tss_ASVs$`16S_ASVs`) %>% as("matrix") %>% as.data.frame(stringsAsFactors = FALSE)
    rownames(tax_table_16S) <- rownames(tax_table(tss_ASVs$`16S_ASVs`))

    tax_table_ITS <- tax_table(tss_ASVs$`ITS_ASVs`) %>% as("matrix") %>% as.data.frame(stringsAsFactors = FALSE)
    rownames(tax_table_ITS) <- rownames(tax_table(tss_ASVs$`ITS_ASVs`))

    # silva
    tax_table(tss_ASVs$`16S_ASVs`) <- tax_table_16S %>% select(starts_with("silva")) %>% as.matrix() %>% tax_table()
    tss_silva <- list(
        `16S_silva_genus` = tss_ASVs$`16S_ASVs` %>% tax_glom("silva_genus"),
        `16S_silva_family` = tss_ASVs$`16S_ASVs` %>% tax_glom("silva_family"),
        `16S_silva_order` = tss_ASVs$`16S_ASVs` %>% tax_glom("silva_order"),
        `16S_silva_class` = tss_ASVs$`16S_ASVs` %>% tax_glom("silva_class"),
        `16S_silva_phylum` = tss_ASVs$`16S_ASVs` %>% tax_glom("silva_phylum")
    )

    # rdp
    tax_table(tss_ASVs$`16S_ASVs`) <- tax_table_16S %>% select(starts_with("rdp")) %>% as.matrix() %>% tax_table()
    tss_rdp <- list(
        `16S_rdp_genus` = tss_ASVs$`16S_ASVs` %>% tax_glom("rdp_genus"),
        `16S_rdp_family` = tss_ASVs$`16S_ASVs` %>% tax_glom("rdp_family"),
        `16S_rdp_order` = tss_ASVs$`16S_ASVs` %>% tax_glom("rdp_order"),
        `16S_rdp_class` = tss_ASVs$`16S_ASVs` %>% tax_glom("rdp_class"),
        `16S_rdp_phylum` = tss_ASVs$`16S_ASVs` %>% tax_glom("rdp_phylum")
    )
    
    # decipher
    tax_table(tss_ASVs$`16S_ASVs`) <- tax_table_16S %>% select(starts_with("decipher")) %>% as.matrix() %>% tax_table()
    tss_decipher <- list(
        `16S_decipher_genus` = tss_ASVs$`16S_ASVs` %>% tax_glom("decipher_genus"),
        `16S_decipher_family` = tss_ASVs$`16S_ASVs` %>% tax_glom("decipher_family"),
        `16S_decipher_order` = tss_ASVs$`16S_ASVs` %>% tax_glom("decipher_order"),
        `16S_decipher_class` = tss_ASVs$`16S_ASVs` %>% tax_glom("decipher_class"),
        `16S_decipher_phylum` = tss_ASVs$`16S_ASVs` %>% tax_glom("decipher_phylum")
    )

    # unite
    tax_table(tss_ASVs$`ITS_ASVs`) <- tax_table_ITS %>% select(starts_with("unite")) %>% as.matrix() %>% tax_table()
    tss_unite <- list(
        `ITS_unite_genus` = tss_ASVs$`ITS_ASVs` %>% tax_glom("unite_genus"),
        `ITS_unite_family` = tss_ASVs$`ITS_ASVs` %>% tax_glom("unite_family"),
        `ITS_unite_order` = tss_ASVs$`ITS_ASVs` %>% tax_glom("unite_order"),
        `ITS_unite_class` = tss_ASVs$`ITS_ASVs` %>% tax_glom("unite_class"),
        `ITS_unite_phylum` = tss_ASVs$`ITS_ASVs` %>% tax_glom("unite_phylum")
    )
    
    # restore tax_tables
    tax_table(tss_ASVs$`16S_ASVs`) <- tax_table_16S %>% as.matrix() %>% tax_table()
    tax_table(tss_ASVs$`ITS_ASVs`) <- tax_table_ITS %>% as.matrix() %>% tax_table()

    # concatenate results
    ps.list.tss <- c(tss_ASVs, tss_silva, tss_rdp, tss_decipher, tss_unite)

# save original (ps.list) and normalized (ps.list.tss) phyloseq objects

    save(ps.list, ps.list.tss, file = filenames$output)
