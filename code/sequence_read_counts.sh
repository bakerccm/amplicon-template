#!/bin/bash

## This needs to be updated. It does not run hands off. ##

OUTPUT_DIR="out/combined/sequence_counts"
rm -rf ${OUTPUT_DIR}
mkdir ${OUTPUT_DIR}

# get read counts at various stages of the pipeline

## prior to demultiplexing

for f in data/*.fastq.gz; do SAMPLE=`basename ${f}`; SEQS=`echo $(zcat ${f} | wc -l)/4|bc`; echo ${SAMPLE} ${SEQS} >>${OUTPUT_DIR}/pre_demultiplexing_seqs.txt; done &

## after demultiplexing

for f in out/16S/demultiplexed/*-R1.fastq.gz; do SAMPLE=`basename ${f} -R1.fastq.gz`; SEQS=`echo $(zcat ${f} | wc -l)/4|bc`; echo ${SAMPLE} ${SEQS} >>${OUTPUT_DIR}/post_demultiplexing_seqs_16S.txt; done &
for f in out/ITS/demultiplexed/*-R1.fastq.gz; do SAMPLE=`basename ${f} -R1.fastq.gz`; SEQS=`echo $(zcat ${f} | wc -l)/4|bc`; echo ${SAMPLE} ${SEQS} >>${OUTPUT_DIR}/post_demultiplexing_seqs_ITS.txt; done &

## after demultiplexing and filtering

for f in out/16S/demultiplexed_no_Ns/*-R1.fastq.gz; do SAMPLE=`basename ${f} -R1.fastq.gz`; SEQS=`echo $(zcat ${f} | wc -l)/4|bc`; echo ${SAMPLE} ${SEQS} >>${OUTPUT_DIR}/post_demultiplexing_noNs_seqs_16S.txt; done &
for f in out/ITS/demultiplexed_no_Ns/*-R1.fastq.gz; do SAMPLE=`basename ${f} -R1.fastq.gz`; SEQS=`echo $(zcat ${f} | wc -l)/4|bc`; echo ${SAMPLE} ${SEQS} >>${OUTPUT_DIR}/post_demultiplexing_noNs_seqs_ITS.txt; done &

## after cutadapt

for f in out/16S/cutadapt/*-R1.fastq.gz; do SAMPLE=`basename ${f} -R1.fastq.gz`; SEQS=`echo $(zcat ${f} | wc -l)/4|bc`; echo ${SAMPLE} ${SEQS} >>${OUTPUT_DIR}/post_cutadapt_seqs_16S.txt; done &
for f in out/ITS/cutadapt/*-R1.fastq.gz; do SAMPLE=`basename ${f} -R1.fastq.gz`; SEQS=`echo $(zcat ${f} | wc -l)/4|bc`; echo ${SAMPLE} ${SEQS} >>${OUTPUT_DIR}/post_cutadapt_seqs_ITS.txt; done &

# post filter and trim

for f in out/16S/filterAndTrim/*-R1.fastq.gz; do SAMPLE=`basename ${f} -R1.fastq.gz`; SEQS=`echo $(zcat ${f} | wc -l)/4|bc`; echo ${SAMPLE} ${SEQS} >>${OUTPUT_DIR}/post_filterAndTrim_seqs_16S.txt; done &
for f in out/ITS/filterAndTrim/*-R1.fastq.gz; do SAMPLE=`basename ${f} -R1.fastq.gz`; SEQS=`echo $(zcat ${f} | wc -l)/4|bc`; echo ${SAMPLE} ${SEQS} >>${OUTPUT_DIR}/post_filterAndTrim_seqs_ITS.txt; done &

# post dada
# use conda activate /mnt/home/ernakovich/cbaker/SERDP-RDX/.snakemake/conda/7e848fb8afaf564ae3c67080436f2ee5_
R
library("dada2")
x <- readRDS("out/16S/dada/sequence_table.rds")
mysums <- apply(x, MAR =1 , sum)
mysums.df <- data.frame(sample = names(mysums), reads = mysums)
rownames(mysums.df) <- NULL
write.table(mysums.df, file = "out/combined/sequence_counts/post_dada_sequencetable_16S.txt", row.names = FALSE, quote= FALSE)

x <- readRDS("out/ITS/dada/sequence_table.rds")
mysums <- apply(x, MAR =1 , sum)
mysums.df <- data.frame(sample = names(mysums), reads = mysums)
rownames(mysums.df) <- NULL
write.table(mysums.df, file = "out/combined/sequence_counts/post_dada_sequencetable_ITS.txt", row.names = FALSE, quote= FALSE)

################################################################################
# post chimeras

R
x <- readRDS("out/16S/remove_chimeras/sequence_table.rds")
mysums <- apply(x, MAR =1 , sum)
mysums.df <- data.frame(sample = names(mysums), reads = mysums)
rownames(mysums.df) <- NULL
write.table(mysums.df, file = "out/combined/sequence_counts/post_remove_chimeras_sequencetable_16S.txt", row.names = FALSE, quote= FALSE)

x <- readRDS("out/ITS/remove_chimeras/sequence_table.rds")
mysums <- apply(x, MAR =1 , sum)
mysums.df <- data.frame(sample = names(mysums), reads = mysums)
rownames(mysums.df) <- NULL
write.table(mysums.df, file = "out/combined/sequence_counts/post_remove_chimeras_sequencetable_ITS.txt", row.names = FALSE, quote= FALSE)

################################################################################
# post export to phyloseq (should be the same)

# use conda activate /mnt/home/ernakovich/cbaker/SERDP-RDX/.snakemake/conda/b15567e953bbde4d6cf579980be0ea53_
R
library(phyloseq)
x <- readRDS("out/16S/phyloseq/phyloseq.rds")
mysums <- sample_sums(x)
mysums.df <- data.frame(sample = names(mysums), reads = mysums)
rownames(mysums.df) <- NULL
write.table(mysums.df, file = "out/combined/sequence_counts/post_export_phyloseq_16S.txt", row.names = FALSE, quote= FALSE)
x <- readRDS("out/16S/phyloseq/phyloseq_cleaned.rds")
mysums <- sample_sums(x)
mysums.df <- data.frame(sample = names(mysums), reads = mysums)
rownames(mysums.df) <- NULL
write.table(mysums.df, file = "out/combined/sequence_counts/post_export_phyloseq_16S_cleaned.txt", row.names = FALSE, quote= FALSE)
x <- readRDS("out/ITS/phyloseq/phyloseq.rds")
mysums <- sample_sums(x)
mysums.df <- data.frame(sample = names(mysums), reads = mysums)
rownames(mysums.df) <- NULL
write.table(mysums.df, file = "out/combined/sequence_counts/post_export_phyloseq_ITS.txt", row.names = FALSE, quote= FALSE)

# post normalize (should be the same)


R
library(phyloseq)
load ("out/combined/serdp_normalized.rdata")
mysums <- sample_sums(ps.list$`16S`)
mysums.df <- data.frame(sample = names(mysums), reads = mysums)
rownames(mysums.df) <- NULL
write.table(mysums.df, file = "out/combined/sequence_counts/post_export_phyloseq_normalized_16S.txt", row.names = FALSE, quote= FALSE)
mysums <- sample_sums(ps.list$`ITS`)
mysums.df <- data.frame(sample = names(mysums), reads = mysums)
rownames(mysums.df) <- NULL
write.table(mysums.df, file = "out/combined/sequence_counts/post_export_phyloseq_normalized_ITS.txt", row.names = FALSE, quote= FALSE)
