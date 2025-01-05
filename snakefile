# snakefile for processing amplicon data

# Chris Baker
# https://github.com/bakerccm
# 1 November 2023

from snakemake.utils import min_version
min_version("7.25.0")

configfile: 'config/config.yaml'

wildcard_constraints:
    sample = '[^/.]+',
    file = '[^/.]+',
    read = 'R[12]'

################################
# get sample metadata

import pandas as pd
import os # for os.path.join()

DATASETS = config['datasets'] # 16S, ITS etc

METADATA = {} # dict of metadata tables
SAMPLES = {} # dict of sample id lists

for dataset in DATASETS:
    METADATA[dataset] = pd.read_csv(config['sample_metadata'][dataset], sep = '\t', index_col = 'SampleID')
    SAMPLES[dataset] = list(METADATA[dataset].index)

################################
# default rule

rule all:
    input:
        expand('out/{dataset}/demultiplexed_qc/multiqc_report.html', dataset = DATASETS),
        [os.path.join('out', dataset, 'read_quality_profiles', sample + '.pdf') 
            for dataset in DATASETS 
            for sample in SAMPLES[dataset]],
        "out/combined/amplicon.rds",
        "out/sequence_counts/summary/aggregated.csv"

################################
# demultiplex samples and run qc on demultiplexed files

rule all_demultiplex_qc:
    input:
        expand('out/{dataset}/demultiplexed_qc/multiqc_report.html', dataset = DATASETS)

# make barcode lists from metadata files
# N.B. pulls out first column (should be sample) and fourth column (should be barcode) but
# does not verify that these columns are correct so use caution if metadata file is reformatted
rule barcode_list:
    input:
        lambda wildcards: config['sample_metadata'][wildcards.dataset]
    output:
        'out/{dataset}/barcodes.txt'
    shell:
        '''
        awk 'BEGIN {{FS="\t";OFS="\t"}} !(NR==1) {{print $1,$4}}' {input} >{output}
        '''

# generates a demultiplex rule for each dataset
# - list of output files varies by dataset but this can't be left as a wildcard in the output expand() if a single rule is used for both datasets
# - alternatively, if you write a rule that just specifies a single sample's outputs rather than all the output files for a dataset, the rule runs once per sample rather than once per dataset
# - see discussion at https://stackoverflow.com/questions/41135801/snakemake-best-practice-for-demultiplexing
# - note: requires snakemake 5.31.0 or later for 'name' keyword to work
for dataset in DATASETS:
    rule:
        name:
            'demultiplex_' + dataset
        input:
            barcodes = os.path.join('out', dataset, 'barcodes.txt'),
            read1 = config['raw_data'][dataset]['read1'],
            read2 = config['raw_data'][dataset]['read2'],
            index = config['raw_data'][dataset]['index']
        output:
            [os.path.join('out', dataset, 'demultiplexed', sample + '-' + read + '.fastq') 
                for sample in SAMPLES[dataset] 
                for read in ['R1','R2']]
        params:
            output_dir = os.path.join('out', dataset, 'demultiplexed')
        conda:
            'envs/illumina-utils-2.12.yaml'
        shell:
            '''
            # delete any existing output files
                for file in {output}
                do
                    rm -f ${{file}}
                done
            # demultiplex raw data files
                iu-demultiplex -s {input.barcodes} --r1 {input.read1} --r2 {input.read2} \
                    --index {input.index} -o {params.output_dir} {config[iu-demultiplex]}
            # create empty files for any samples that did not have reads
                for file in {output}
                do
                    if [ ! -f ${{file}} ]
                    then
                        touch ${{file}}
                    fi
                done
            '''

rule compress_demultiplexed_fastq:
    input:
        'out/{dataset}/demultiplexed/{sample}-{read}.fastq'
    output:
        'out/{dataset}/demultiplexed/{sample}-{read}.fastq.gz'
    shell:
        'gzip {input}'

rule demultiplex_fastqc:
    input:
        'out/{dataset}/demultiplexed/{sample}-{read}.fastq.gz'
    output:
        'out/{dataset}/demultiplexed_qc/{sample}-{read}_fastqc.html'
    params:
        output_dir='out/{dataset}/demultiplexed_qc'
    conda:
        'envs/fastqc-0.11.9.yaml'
    shell:
        'fastqc -o {params.output_dir} {input}'

# N.B. this rule asks for fastqc files for all demultiplexed samples
rule demultiplex_multiqc:
    input:
        fastqc = lambda wildcards: [os.path.join('out', wildcards.dataset, 'demultiplexed_qc', sample + '-' + read + '_fastqc.html') 
            for sample in SAMPLES[wildcards.dataset]
            for read in ['R1','R2']]
    output:
        'out/{dataset}/demultiplexed_qc/multiqc_report.html'
    params:
        inputdir = 'out/{dataset}/demultiplexed_qc',
        outputdir = 'out/{dataset}/demultiplexed_qc'
    log:
        'out/{dataset}/demultiplexed_qc/multiqc_report.log'
    conda:
        'envs/multiqc-1.22.yaml'
    shell:
        'multiqc {config[multiqc]} --force -o {params.outputdir} {params.inputdir} 2>{log}'

################################
# remove Ns and check primer orientation in preparation for running cutadapt

rule all_filter_Ns:
    input:
        [os.path.join('out', dataset, 'demultiplexed_no_Ns', sample + '-' + file) 
            for dataset in DATASETS 
            for sample in SAMPLES[dataset] 
            for file in ['R1.fastq.gz','R2.fastq.gz','summary.txt']]

rule filter_Ns:
    input:
        read1='out/{dataset}/demultiplexed/{sample}-R1.fastq.gz',
        read2='out/{dataset}/demultiplexed/{sample}-R2.fastq.gz'
    output:
        read1='out/{dataset}/demultiplexed_no_Ns/{sample}-R1.fastq.gz',
        read2='out/{dataset}/demultiplexed_no_Ns/{sample}-R2.fastq.gz',
        summary='out/{dataset}/demultiplexed_no_Ns/{sample}-summary.txt'
    params:
        primer_fwd=lambda wildcards: config['primers'][wildcards.dataset]['forward'],
        primer_rev=lambda wildcards: config['primers'][wildcards.dataset]['reverse']
    conda:
        'envs/dada2-1.18.0.yaml'
    shell:
        '''
        # test input file and skip processing if empty
        if [ $(zcat {input.read1} | wc -l) -eq 0 ]
        then
            OUTPUT_READ1={output.read1} &&
                touch ${{OUTPUT_READ1%.gz}} &&
                gzip ${{OUTPUT_READ1%.gz}}
            OUTPUT_READ2={output.read2} &&
                touch ${{OUTPUT_READ2%.gz}} &&
                gzip ${{OUTPUT_READ2%.gz}}
            touch {output.summary}
        else
            Rscript code/filter_Ns_and_check_primers.R \
                {input.read1} {input.read2} \
                {params.primer_fwd} {params.primer_rev} \
                {output.read1} {output.read2} {output.summary}
        fi                                                                                           
        '''


################################
# use cutadapt to remove primers
# note don't length trim ITS sequences

rule all_cutadapt:
    input:
        [os.path.join('out', dataset, 'cutadapt', sample + '-' + file) 
            for dataset in DATASETS 
            for sample in SAMPLES[dataset] 
            for file in ['R1.fastq.gz','R2.fastq.gz']]

rule cutadapt:
    input:
        read1='out/{dataset}/demultiplexed_no_Ns/{sample}-R1.fastq.gz',
        read2='out/{dataset}/demultiplexed_no_Ns/{sample}-R2.fastq.gz'
    output:
        read1='out/{dataset}/cutadapt/{sample}-R1.fastq.gz',
        read2='out/{dataset}/cutadapt/{sample}-R2.fastq.gz'
    log:
        'out/{dataset}/cutadapt/{sample}.log'
    params:
        primer_fwd = lambda wildcards: config['primers'][wildcards.dataset]['forward'],
        primer_rev = lambda wildcards: config['primers'][wildcards.dataset]['reverse'],
        min_length = config['cutadapt']['min_length']
    conda:
        'envs/cutadapt-4.6.yaml'
    shell:
        '''
        # test input file and skip processing if empty
        if [ $(zcat {input.read1} | wc -l) -eq 0 ]
        then
            OUTPUT_READ1={output.read1} &&
                touch ${{OUTPUT_READ1%.gz}} &&
                gzip ${{OUTPUT_READ1%.gz}}
            OUTPUT_READ2={output.read2} &&
                touch ${{OUTPUT_READ2%.gz}} &&
                gzip ${{OUTPUT_READ2%.gz}}
        else
            # these are the primer sequences
                primer_fwd={params.primer_fwd}
                primer_rev={params.primer_rev}
            # get reverse complements
                primer_fwd_rc=`echo ${{primer_fwd}} | tr ACGTMRWSYKVHDBNacgtmrwsykvhdbn TGCAKYWSRMBDHVNtgcakywsrmbdhvn | rev`
                primer_rev_rc=`echo ${{primer_rev}} | tr ACGTMRWSYKVHDBNacgtmrwsykvhdbn TGCAKYWSRMBDHVNtgcakywsrmbdhvn | rev`
            cutadapt -g ${{primer_fwd}} -a ${{primer_rev_rc}} -G ${{primer_rev}} -A ${{primer_fwd_rc}} -n 2 -m {params.min_length} \
                -o {output.read1} -p {output.read2} {input.read1} {input.read2} >{log}
        fi
        '''

################################
# plot and inspect read quality profiles

rule all_read_quality_profiles:
    input:
        [os.path.join('out', dataset, 'read_quality_profiles', sample + '.pdf') 
            for dataset in DATASETS 
            for sample in SAMPLES[dataset]]

rule read_quality_profiles:
    input:
        read1='out/{dataset}/cutadapt/{sample}-R1.fastq.gz',
        read2='out/{dataset}/cutadapt/{sample}-R2.fastq.gz'
    output:
        'out/{dataset}/read_quality_profiles/{sample}.pdf'
    conda:
        'envs/dada2-1.18.0.yaml'
    shell:
        'Rscript code/plot_read_profile.R {input.read1} {input.read2} {output}'

################################
# consider using Figaro to assess trimming lengths
# https://github.com/Zymo-Research/figaro#figaro

################################
# filter and trim (in dada2)

# note use of different rules for 16S and ITS
# - 16S workflow uses truncLen to force fixed read lengths
# - ITS workflow does not do this (since ITS can vary in length), but does enforce minLen to remove spurious short reads

rule all_filter_and_trim:
    input:
        [os.path.join('out', dataset, 'filterAndTrim', sample + '-' + read + '.fastq.gz') 
            for dataset in DATASETS 
            for sample in SAMPLES[dataset] 
            for read in ['R1','R2']]

rule filter_and_trim_16S:
    input:
        read1 = 'out/16S/cutadapt/{sample}-R1.fastq.gz',
        read2 = 'out/16S/cutadapt/{sample}-R2.fastq.gz'
    output:
        read1 = 'out/16S/filterAndTrim/{sample}-R1.fastq.gz',
        read2 = 'out/16S/filterAndTrim/{sample}-R2.fastq.gz'
    params:
        maxN = config['fastqPairedFilter']['16S']['maxN'],
        truncQ = config['fastqPairedFilter']['16S']['truncQ'],
        maxEE_read1 = config['fastqPairedFilter']['16S']['maxEE_read1'],
        maxEE_read2 = config['fastqPairedFilter']['16S']['maxEE_read2'],
        truncLen_read1 = config['fastqPairedFilter']['16S']['truncLen_read1'],
        truncLen_read2 = config['fastqPairedFilter']['16S']['truncLen_read2']
    conda:
        'envs/dada2-1.18.0.yaml'
    shell:
        '''
        # test input file and skip processing if empty
        if [ $(zcat {input.read1} | wc -l) -eq 0 ]
        then
            OUTPUT_READ1={output.read1} &&
                touch ${{OUTPUT_READ1%.gz}} &&
                gzip ${{OUTPUT_READ1%.gz}}
            OUTPUT_READ2={output.read2} &&
                touch ${{OUTPUT_READ2%.gz}} &&
                gzip ${{OUTPUT_READ2%.gz}}
        else
            Rscript code/filter_and_trim_16S.R \
                {input.read1} {input.read2} {output.read1} {output.read2} \
                {params.maxN} {params.truncQ} {params.maxEE_read1} {params.maxEE_read2} {params.truncLen_read1} {params.truncLen_read2}
        fi
        '''

rule filter_and_trim_ITS:
    input:
        read1 = 'out/ITS/cutadapt/{sample}-R1.fastq.gz',
        read2 = 'out/ITS/cutadapt/{sample}-R2.fastq.gz'
    output:
        read1 = 'out/ITS/filterAndTrim/{sample}-R1.fastq.gz',
        read2 = 'out/ITS/filterAndTrim/{sample}-R2.fastq.gz'
    params:
        maxN = config['fastqPairedFilter']['ITS']['maxN'],
        truncQ = config['fastqPairedFilter']['ITS']['truncQ'],
        maxEE_read1 = config['fastqPairedFilter']['ITS']['maxEE_read1'],
        maxEE_read2 = config['fastqPairedFilter']['ITS']['maxEE_read2'],
        minLen = config['fastqPairedFilter']['ITS']['minLen']
    conda:
        'envs/dada2-1.18.0.yaml'
    shell:
        '''
        # test input file and skip processing if empty
        if [ $(zcat {input.read1} | wc -l) -eq 0 ]
        then
            OUTPUT_READ1={output.read1} &&
                touch ${{OUTPUT_READ1%.gz}} &&
                gzip ${{OUTPUT_READ1%.gz}}
            OUTPUT_READ2={output.read2} &&
                touch ${{OUTPUT_READ2%.gz}} &&
                gzip ${{OUTPUT_READ2%.gz}}
        else
            Rscript code/filter_and_trim_ITS.R \
                {input.read1} {input.read2} {output.read1} {output.read2} \
                {params.maxN} {params.truncQ} {params.maxEE_read1} {params.maxEE_read2} {params.minLen}
        fi
        '''

################################
# learn error rates in dada2

rule all_learnerrors_inputs:
    input:
        [os.path.join('out', dataset, 'learnerrors', read, sample + '-' + read + '.fastq.gz') 
            for dataset in DATASETS 
            for sample in SAMPLES[dataset] 
            for read in ['R1','R2']]

rule learnerrors_inputs:
    input:
        read1 = 'out/{dataset}/filterAndTrim/{sample}-R1.fastq.gz',
        read2 = 'out/{dataset}/filterAndTrim/{sample}-R2.fastq.gz'
    output:
        read1 = 'out/{dataset}/learnerrors/R1/{sample}-R1.fastq.gz',
        read2 = 'out/{dataset}/learnerrors/R2/{sample}-R2.fastq.gz'
    shell:
        '''
        ln -s ../../../../{input.read1} {output.read1}
        ln -s ../../../../{input.read2} {output.read2}
        '''

rule all_learnerrors:
    input:
        # one job for each read direction in each dataset
        [os.path.join('out', dataset, 'learnerrors', read + '_error_profile.rds') 
            for dataset in DATASETS 
            for read in ['R1','R2']]

rule learnerrors:
    input:
        # all the fastq files for the specified read direction (R1 or R2) for the specified dataset (16S or ITS)
        lambda wildcards: [os.path.join('out', wildcards.dataset, 'learnerrors', wildcards.read, sample + '-' + wildcards.read + '.fastq.gz') 
            for sample in SAMPLES[wildcards.dataset]]
    output:
        rds = 'out/{dataset}/learnerrors/{read}_error_profile.rds',
        pdf = 'out/{dataset}/learnerrors/{read}_error_profile.pdf'
    log:
        'out/{dataset}/learnerrors/{read}_error_profile.log'
    params:
        input_dir = 'out/{dataset}/learnerrors/{read}',
        nbases = config['learnErrors']['nbases']
    conda:
        'envs/dada2-1.18.0.yaml'
    threads: 1
    shell:
        'Rscript code/learn_errors.R {params.input_dir} {output.rds} {output.pdf} {params.nbases} 1>{log}'

################################
# sequence inference for each sample
# (dereplicate --> run dada --> merge pairs)

# consider running dada in pooled mode to pick up rare sequences
# see https://benjjneb.github.io/dada2/tutorial.html

rule all_dada:
    input:
        [os.path.join('out', dataset, 'dada', sample + '.rds') 
            for dataset in DATASETS 
            for sample in SAMPLES[dataset]]

# this rule appears to be memory intensive, especially for samples with many reads e.g the respiration samples
# -- it may be necessary to limit execution to a couple of jobs at a time to ensure enough memory is available
rule dada:
    input:
        fastq_read1 = 'out/{dataset}/learnerrors/R1/{sample}-R1.fastq.gz',
        fastq_read2 = 'out/{dataset}/learnerrors/R2/{sample}-R2.fastq.gz',
        error_profile_read1 = 'out/{dataset}/learnerrors/R1_error_profile.rds',
        error_profile_read2 = 'out/{dataset}/learnerrors/R2_error_profile.rds'
    output:
        'out/{dataset}/dada/{sample}.rds'
    log:
        'out/{dataset}/dada/{sample}.log'
    conda:
        'envs/dada2-1.18.0.yaml'
    shell:
        '''
        # test input file and skip processing if empty
        if [ $(zcat {input.fastq_read1} | wc -l) -eq 0 ]
        then
            touch {output}
        else
            Rscript code/dada.R {input.fastq_read1} {input.fastq_read2} \
                {input.error_profile_read1} {input.error_profile_read2} {output} >{log}
        fi
        '''

rule all_dada_merge:
    input:
        expand('out/{dataset}/dada-merge/sequence_table.rds', dataset = DATASETS)

# each dataset runs in ~1min with a single core
rule dada_merge_samples:
    input:
        read1_fastq_files = lambda wildcards: [os.path.join('out', wildcards.dataset, 'learnerrors', 'R1', sample + '-R1.fastq.gz') 
            for sample in SAMPLES[wildcards.dataset]],
        rds_files = lambda wildcards: [os.path.join('out', wildcards.dataset, 'dada', sample + '.rds') 
            for sample in SAMPLES[wildcards.dataset]]
    output:
        sample_list = 'out/{dataset}/dada-merge/sample_list.txt',
        sequence_table = 'out/{dataset}/dada-merge/sequence_table.rds'
    log:
        'out/{dataset}/dada-merge/sequence_table.log'
    conda:
        'envs/dada2-1.18.0.yaml'
    shell:
        '''
        # write list of non-empty input files to file
            rm -rf {output.sample_list}
            for READ1_FASTQ_FILE in {input.read1_fastq_files}
            do
                SAMPLE=`basename ${{READ1_FASTQ_FILE%-R1.fastq.gz}}`
                RDS_FILE="out/{wildcards.dataset}/dada/${{SAMPLE}}.rds"
                if [ $(zcat ${{READ1_FASTQ_FILE}} | wc -l) -gt 0 ]
                then
                    echo ${{RDS_FILE}} >>{output.sample_list}
                fi
            done
        # make sequence table from non-empty samples
            Rscript code/merge_samples.R {output.sample_list} {output.sequence_table} >{log}
        '''

################################
# merge runs (see https://benjjneb.github.io/dada2/bigdata_paired.html)
# -- if applicable

################################
# remove chimeras

rule all_remove_chimeras:
    input:
        expand('out/{dataset}/remove_chimeras/sequence_table.rds', dataset = DATASETS)

# each dataset runs in ~1min with a single core
rule remove_chimeras:
    input:
        'out/{dataset}/dada-merge/sequence_table.rds'
    output:
        'out/{dataset}/remove_chimeras/sequence_table.rds'
    log:
        'out/{dataset}/remove_chimeras/sequence_table.log'
    conda:
        'envs/dada2-1.18.0.yaml'
    shell:
        'Rscript code/remove_chimeras.R {input} {output} >{log}'

################################
# assign taxonomy

rule all_taxonomy:
    input:
        expand('out/16S/assign_taxonomy/{database}_{matches}.rds', database = ['rdp', 'silva'], matches = ['single', 'multiple']),
        'out/ITS/assign_taxonomy/unite.rds', 'out/16S/decipher/silva.rds'

# uses dada2::assignTaxonomy() and dada2::assignSpecies()
# returns a file of unambiguous matches and a file of multiple matches
# - takes ~15 min to run RDP and ~25 min to run Silva with 2 cores, 16GB per core
rule assign_taxonomy_16S:
    input:
        sequence_table = 'out/16S/remove_chimeras/sequence_table.rds',
        assignTaxonomy_refFasta = lambda wildcards: config['databases']['assignTaxonomy'][wildcards.database],
        assignSpecies_refFasta = lambda wildcards: config['databases']['assignSpecies'][wildcards.database]
    output:
        single_match = 'out/16S/assign_taxonomy/{database}_single.rds',
        multiple_matches = 'out/16S/assign_taxonomy/{database}_multiple.rds'
    log:
        'out/16S/assign_taxonomy/{database}.log'
    conda:
        'envs/dada2-1.18.0.yaml'
    threads: config['assignTaxonomy']['threads']
    shell:
        '''
        Rscript code/assign_taxonomy_16S.R {input.sequence_table} \
            {input.assignTaxonomy_refFasta} {input.assignSpecies_refFasta} \
            {output.single_match} {output.multiple_matches} >{log}
        '''

# uses dada2::assignTaxonomy() only
rule assign_taxonomy_ITS:
    input:
        sequence_table = 'out/ITS/remove_chimeras/sequence_table.rds',
        assignTaxonomy_refFasta = lambda wildcards: config['databases']['assignTaxonomy'][wildcards.database]
    output:
        'out/ITS/assign_taxonomy/{database}.rds'
    log:
        'out/ITS/assign_taxonomy/{database}.log'
    conda:
        'envs/dada2-1.18.0.yaml'
    threads: config['assignTaxonomy']['threads']
    shell:
        'Rscript code/assign_taxonomy_ITS.R {input.sequence_table} {input.assignTaxonomy_refFasta} {output} >{log}'

# uses DECIPHER::IdTaxa()
rule decipher:
    input:
        sequence_table = 'out/16S/remove_chimeras/sequence_table.rds',
        decipher_ref = lambda wildcards: config['databases']['DECIPHER'][wildcards.database]
    output:
        'out/16S/decipher/{database}.rds'
    log:
        'out/16S/decipher/{database}.log'
    conda:
        'envs/dada2-1.18.0-decipher.yaml'
    threads: config['decipher']['threads']
    shell:
        'Rscript code/decipher.R {input.sequence_table} {input.decipher_ref} {output} >{log}'

################################
# export as phyloseq object

rule all_export_to_phyloseq:
    input:
        expand("out/{dataset}/phyloseq/phyloseq.rds", dataset = DATASETS)

rule export_to_phyloseq_16S:
    input:
        sequence_table = "out/16S/remove_chimeras/sequence_table.rds",
        silva_single = "out/16S/assign_taxonomy/silva_single.rds",
        silva_multiple = "out/16S/assign_taxonomy/silva_multiple.rds",
        rdp_single = "out/16S/assign_taxonomy/rdp_single.rds",
        rdp_multiple = "out/16S/assign_taxonomy/rdp_multiple.rds",
        decipher = "out/16S/decipher/silva.rds",
        metadata = config['sample_metadata']['16S']
    output:
        phyloseq="out/16S/phyloseq/phyloseq.rds",
        fasta="out/16S/phyloseq/phyloseq.fasta"
    log:
        "out/16S/phyloseq/phyloseq.log"
    conda:
        "envs/phyloseq.yaml"
    shell:
        '''
        Rscript code/export_to_phyloseq_16S.R {input.sequence_table} \
        {input.silva_single} {input.silva_multiple} {input.rdp_single} {input.rdp_multiple} {input.decipher} \
        {input.metadata} {output.phyloseq} {output.fasta} >{log}
        '''

rule export_to_phyloseq_ITS:
    input:
        sequence_table = "out/ITS/remove_chimeras/sequence_table.rds",
        unite = "out/ITS/assign_taxonomy/unite.rds",
        metadata = config['sample_metadata']['ITS']
    output:
        "out/ITS/phyloseq/phyloseq.rds"
    log:
        "out/ITS/phyloseq/phyloseq.log"
    conda:
        "envs/phyloseq.yaml"
    shell:
        "Rscript code/export_to_phyloseq_ITS.R {input.sequence_table} {input.unite} {input.metadata} {output} >{log}"

################################
rule remove_chloroplasts:
    input:
        rds = "out/16S/phyloseq/phyloseq.rds"
    output:
        rds = "out/16S/phyloseq/phyloseq_cleaned.rds",
        fasta = "out/16S/phyloseq/phyloseq_cleaned.fasta"
    conda:
        "envs/phyloseq-tidyverse.yaml"
    shell:
        "Rscript code/remove_chloroplasts.R {input.rds} {output.rds} {output.fasta}"

################################
# Note: SEPP needs to be installed separately (not through conda) in the `/software/sepp-package` folder.
# The `sepp_env` conda environment provides a suitable runtime environment but does not actually install SEPP.
# The script `scripts/install_sepp.sh` provides suitable commands for installing SEPP.
rule run_sepp:
    input:
        "out/16S/phyloseq/phyloseq_cleaned.fasta"
    output:
        "out/16S/phyloseq/phyloseq_cleaned_placement.tog.tre"
    conda:
        "envs/sepp_env.yaml"
    threads:
        config['sepp']['threads']
    shell:
        # running sepp by cd'ing into this directory seems to be the way to do it but it requires some convoluted relative paths to get it to work here
        '''
        cd out/16S/phyloseq
        ../../../software/sepp-package/run-sepp.sh phyloseq_cleaned.fasta phyloseq_cleaned -x {threads}
        '''

################################
rule combine_data:
    input:
        phyloseq_16S = "out/16S/phyloseq/phyloseq_cleaned.rds",
        sepp_tree_16S = "out/16S/phyloseq/phyloseq_cleaned_placement.tog.tre",
        phyloseq_ITS = "out/ITS/phyloseq/phyloseq.rds"
    output:
        "out/combined/amplicon.rds"
    log:
        "out/combined/amplicon.log"
    conda:
       	"envs/phyloseq-tidyverse.yaml"
    shell:
        "Rscript code/combine_data.R {input.phyloseq_16S} {input.sepp_tree_16S} {input.phyloseq_ITS} {output} >{log}"

################################
# get sequence counts at various stages of the pipeline

# sequence counts prior to demultiplexing

rule sequence_counts_before_demultiplexing:
    input:
        lambda wildcards: config['raw_data'][wildcards.dataset]['read1']
    output:
        "out/sequence_counts/before_demultiplexing/{dataset}/{dataset}.csv"
    shell:
        '''
        DATASET={wildcards.dataset}
        STAGE='before demultiplexing'
        SAMPLE=`basename {input}`
        SEQS=`echo $(zcat {input} | wc -l)/4|bc`
        echo ${{DATASET}},${{STAGE}},${{SAMPLE}},${{SEQS}} >{output}
        '''

rule aggregate_sequence_counts_before_demultiplexing:
    input:
        lambda wildcards: [os.path.join("out", "sequence_counts", "before_demultiplexing", dataset, dataset + ".csv") for dataset in DATASETS]
    output:
        "out/sequence_counts/before_demultiplexing/aggregated.csv"
    shell:
        '''
        echo "dataset,stage,sample,reads" >{output}
        for f in {input}; do cat $f >>{output}; done
        '''

# sequence counts after demultiplexing

rule sequence_counts_after_demultiplexing:
    input:
        "out/{dataset}/demultiplexed/{sample}-R1.fastq.gz"
    output:
        "out/sequence_counts/after_demultiplexing/{dataset}/{sample}.csv"
    shell:
        '''
        DATASET={wildcards.dataset}
        STAGE='after demultiplexing'
        SAMPLE=`basename {input} -R1.fastq.gz`
        SEQS=`echo $(zcat {input} | wc -l)/4|bc`
        echo ${{DATASET}},${{STAGE}},${{SAMPLE}},${{SEQS}} >{output}
        '''

rule aggregate_sequence_counts_after_demultiplexing:
    input:
        lambda wildcards: [os.path.join("out", "sequence_counts", "after_demultiplexing", dataset, sample + ".csv") for dataset in DATASETS for sample in SAMPLES[dataset]]
    output:
        "out/sequence_counts/after_demultiplexing/aggregated.csv"
    shell:
        '''
        echo "dataset,stage,sample,reads" >{output}
        for f in {input}; do cat $f >>{output}; done
        '''

# sequence counts after filtering Ns

rule sequence_counts_after_filtering_Ns:
    input:
        "out/{dataset}/demultiplexed_no_Ns/{sample}-R1.fastq.gz"
    output:
        "out/sequence_counts/after_filtering_Ns/{dataset}/{sample}.csv"
    shell:
        '''
        DATASET={wildcards.dataset}
        STAGE='after filtering Ns'
        SAMPLE=`basename {input} -R1.fastq.gz`
        SEQS=`echo $(zcat {input} | wc -l)/4|bc`
        echo ${{DATASET}},${{STAGE}},${{SAMPLE}},${{SEQS}} >{output}
        '''

rule aggregate_sequence_counts_after_filtering_Ns:
    input:
        lambda wildcards: [os.path.join("out", "sequence_counts", "after_filtering_Ns", dataset, sample + ".csv") for dataset in DATASETS for sample in SAMPLES[dataset]]
    output:
        "out/sequence_counts/after_filtering_Ns/aggregated.csv"
    shell:
        '''
        echo "dataset,stage,sample,reads" >{output}
        for f in {input}; do cat $f >>{output}; done
        '''

# sequence counts after running cutadapt

rule sequence_counts_after_cutadapt:
    input:
        "out/{dataset}/cutadapt/{sample}-R1.fastq.gz"
    output:
        "out/sequence_counts/after_cutadapt/{dataset}/{sample}.csv"
    shell:
        '''
        DATASET={wildcards.dataset}
        STAGE='after cutadapt'
        SAMPLE=`basename {input} -R1.fastq.gz`
        SEQS=`echo $(zcat {input} | wc -l)/4|bc`
        echo ${{DATASET}},${{STAGE}},${{SAMPLE}},${{SEQS}} >{output}
        '''

rule aggregate_sequence_counts_after_cutadapt:
    input:
        lambda wildcards: [os.path.join("out", "sequence_counts", "after_cutadapt", dataset, sample + ".csv") for dataset in DATASETS for sample in SAMPLES[dataset]]
    output:
        "out/sequence_counts/after_cutadapt/aggregated.csv"
    shell:
        '''
        echo "dataset,stage,sample,reads" >{output}
        for f in {input}; do cat $f >>{output}; done
        '''

# sequence counts after filterAndTrim

rule sequence_counts_after_filterAndTrim:
    input:
        "out/{dataset}/filterAndTrim/{sample}-R1.fastq.gz"
    output:
        "out/sequence_counts/after_filterAndTrim/{dataset}/{sample}.csv"
    shell:
        '''
        DATASET={wildcards.dataset}
        STAGE='after filterAndTrim'
        SAMPLE=`basename {input} -R1.fastq.gz`
        SEQS=`echo $(zcat {input} | wc -l)/4|bc`
        echo ${{DATASET}},${{STAGE}},${{SAMPLE}},${{SEQS}} >{output}
        '''

rule aggregate_sequence_counts_after_filterAndTrim:
    input:
        lambda wildcards: [os.path.join("out", "sequence_counts", "after_filterAndTrim", dataset, sample + ".csv") for dataset in DATASETS for sample in SAMPLES[dataset]]
    output:
        "out/sequence_counts/after_filterAndTrim/aggregated.csv"
    shell:
        '''
        echo "dataset,stage,sample,reads" >{output}
        for f in {input}; do cat $f >>{output}; done
        '''

# sequence counts after running dada

rule sequence_counts_after_dada:
    input:
        "out/{dataset}/dada-merge/sequence_table.rds"
    output:
        "out/sequence_counts/after_dada/{dataset}/{dataset}.csv"
    conda:
        "envs/dada2-1.18.0.yaml"
    shell:
        "Rscript code/sequence_counts_after_dada.R {input} {output} {wildcards.dataset} 'after dada'"

rule aggregate_sequence_counts_after_dada:
    input:
        expand("out/sequence_counts/after_dada/{dataset}/{dataset}.csv", dataset = DATASETS)
    output:
        "out/sequence_counts/after_dada/aggregated.csv"
    shell:
        '''
        echo "dataset,stage,sample,reads" >{output}
        for f in {input}; do tail -n +2 $f >>{output}; done
        '''

# sequence counts after removing chimeras (same Rscript as after_dada)

rule sequence_counts_after_removing_chimeras:
    input:
        "out/{dataset}/remove_chimeras/sequence_table.rds"
    output:
        "out/sequence_counts/after_removing_chimeras/{dataset}/{dataset}.csv"
    conda:
        "envs/dada2-1.18.0.yaml"
    shell:
        "Rscript code/sequence_counts_after_dada.R {input} {output} {wildcards.dataset} 'after removing chimeras'"

rule aggregate_sequence_counts_after_removing_chimeras:
    input:
        expand("out/sequence_counts/after_removing_chimeras/{dataset}/{dataset}.csv", dataset = DATASETS)
    output:
        "out/sequence_counts/after_removing_chimeras/aggregated.csv"
    shell:
        '''
        echo "dataset,stage,sample,reads" >{output}
        for f in {input}; do tail -n +2 $f >>{output}; done
        '''

# sequence counts after cleaning

rule sequence_counts_after_cleaning:
    input:
        lambda wildcards: os.path.join("out", wildcards.dataset, "phyloseq", "phyloseq_cleaned.rds") 
            if wildcards.dataset == "16S" 
            else os.path.join("out", wildcards.dataset, "phyloseq", "phyloseq.rds")
    output:
        "out/sequence_counts/after_cleaning/{dataset}/{dataset}.csv"
    conda:
        "envs/phyloseq.yaml"
    shell:
        "Rscript code/sequence_counts_after_cleaning.R {input} {output} {wildcards.dataset} 'after cleaning'"

rule aggregate_sequence_counts_after_cleaning:
    input:
        expand("out/sequence_counts/after_cleaning/{dataset}/{dataset}.csv", dataset = DATASETS)
    output:
        "out/sequence_counts/after_cleaning/aggregated.csv"
    shell:
        '''
        echo "dataset,stage,sample,reads" >{output}
        for f in {input}; do tail -n +2 $f >>{output}; done
        '''

# sequence counts after combining (should be the same as post-cleaning)

rule sequence_counts_after_combining:
    input:
        "out/combined/amplicon.rds"
    output:
        "out/sequence_counts/after_combining/aggregated.csv"
    conda:
        "envs/phyloseq.yaml"
    shell:
        "Rscript code/sequence_counts_after_combining.R {input} {output} 'after combining'"

# collate sequence counts

rule aggregate_sequence_counts_all_stages:
    input:
        before_demultiplexing = "out/sequence_counts/before_demultiplexing/aggregated.csv",
        after_demultiplexing = "out/sequence_counts/after_demultiplexing/aggregated.csv",
        after_filtering_Ns = "out/sequence_counts/after_filtering_Ns/aggregated.csv",
        after_cutadapt = "out/sequence_counts/after_cutadapt/aggregated.csv",
        after_filterAndTrim = "out/sequence_counts/after_filterAndTrim/aggregated.csv",
        after_dada = "out/sequence_counts/after_dada/aggregated.csv",
        after_removing_chimeras = "out/sequence_counts/after_removing_chimeras/aggregated.csv",
        after_cleaning = "out/sequence_counts/after_cleaning/aggregated.csv",
        after_combining = "out/sequence_counts/after_combining/aggregated.csv"
    output:
        "out/sequence_counts/summary/aggregated.csv"
    conda:
        # this rule does not actually use phyloseq but this saves creating a separate environment with here and tidyverse
        "envs/phyloseq-tidyverse.yaml"
    shell:
        '''
        echo "dataset,stage,sample,reads" >{output}
        for f in {input}; do tail -n +2 $f >>{output}; done
        '''

################################
