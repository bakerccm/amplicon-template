# snakefile for processing amplicon data

# Chris Baker
# https://github.com/bakerccm
# 7 April 2025

from snakemake.utils import min_version
min_version('8.27.1')

configfile: 'config/config.yaml'

wildcard_constraints:
    sample = '[^/.]+',
    file = '[^/.]+',
    read = 'R[12]'

################################
# get sample metadata

import pandas as pd
import os # for os.path.join()

# get input file paths from config file
# structure of FILEPATHS should be: marker | file
FILEPATHS = config['filepaths'] # nested dict containing file paths to raw data and metadata files

# get markers from FILEPATHS
MARKERS = list(FILEPATHS.keys())

# get sample metadata tables (METADATA) and sample names (SAMPLES) from files
METADATA = {}
SAMPLES = {}
for marker in MARKERS:
    METADATA[marker] = pd.read_csv(FILEPATHS[marker]['metadata'], sep = '\t', index_col = 'SampleID')
    SAMPLES[marker] = list(METADATA[marker].index)

################################
# default rule

rule all:
    input:
        expand('out/demultiplex_qc/{marker}/multiqc_report.html', marker = MARKERS),
        [os.path.join('out', 'read_quality_profiles', marker, sample + '.pdf') 
            for marker in MARKERS 
            for sample in SAMPLES[marker]],
        'out/phyloseq/phyloseq_list.rds',
        'out/sequence_counts/aggregated.csv'

################################
# demultiplex samples and run qc on demultiplexed files

rule all_demultiplex_qc:
    input:
        expand('out/demultiplex_qc/{marker}/multiqc_report.html', marker = MARKERS)

# make barcode lists from metadata files
# N.B. pulls out sample ID and barcode columns based on headers in first line
# -- ensure that labels for these columns under params are correct
rule barcode_list:
    input:
        lambda wildcards: FILEPATHS[wildcards.marker]['metadata']
    output:
        'out/demultiplex/{marker}/barcodes.txt'
    params:
        sample_id_header = 'SampleID',
        barcode_header = 'BarcodeSequence'
    shell:
        '''
        sam={params.sample_id_header}
        bar={params.barcode_header}
        awk -v sam=${{sam}} -v bar=${{bar}} '
            BEGIN {{FS="\t";OFS="\t"}}
            # get column numbers for sample ID and barcode from header line
            NR==1 {{
                for (i=1; i<=NF; i++) {{
                    if ($i==sam) {{samnum=i}}
                    else if ($i==bar) {{barnum=i}}
                }}
            }}
            # print sample ID and barcode columns from remaining lines of file
            NR>1 {{print $samnum,$barnum}} 
            ' {input} >{output}
        '''

# generates a demultiplex rule for each marker
# - list of output files varies by marker but this can't be left as a wildcard in the output expand() if a single rule is used for both markers
# - alternatively, if you write a rule that just specifies a single sample's outputs rather than all the output files for a marker, the rule runs once per sample rather than once per marker
# - see discussion at https://stackoverflow.com/questions/41135801/snakemake-best-practice-for-demultiplexing
# - note: requires snakemake 5.31.0 or later for 'name' keyword to work
for marker in MARKERS:
    rule:
        name:
            'demultiplex_' + marker
        input:
            barcodes = os.path.join('out', 'demultiplex', marker, 'barcodes.txt'),
            read1 = FILEPATHS[marker]['read1'],
            read2 = FILEPATHS[marker]['read2'],
            index = FILEPATHS[marker]['index']
        output:
            [os.path.join('out', 'demultiplex', marker, sample + '-' + read + '.fastq') 
                for sample in SAMPLES[marker] 
                for read in ['R1','R2']]
        params:
            output_dir = os.path.join('out', 'demultiplex', marker)
        conda:
            'envs/illumina-utils-2.13.yaml'
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
        'out/demultiplex/{marker}/{sample}-{read}.fastq'
    output:
        'out/demultiplex/{marker}/{sample}-{read}.fastq.gz'
    shell:
        'gzip {input}'

rule demultiplex_fastqc:
    input:
        'out/demultiplex/{marker}/{sample}-{read}.fastq.gz'
    output:
        'out/demultiplex_qc/{marker}/{sample}-{read}_fastqc.html'
    params:
        output_dir='out/demultiplex_qc/{marker}'
    conda:
        'envs/fastqc-0.12.1.yaml'
    shell:
        'fastqc -o {params.output_dir} {input}'

# N.B. this rule asks for fastqc files for all demultiplexed samples
rule demultiplex_multiqc:
    input:
        fastqc = lambda wildcards: [os.path.join('out', 'demultiplex_qc', wildcards.marker, sample + '-' + read + '_fastqc.html') 
            for sample in SAMPLES[wildcards.marker]
            for read in ['R1','R2']]
    output:
        'out/demultiplex_qc/{marker}/multiqc_report.html'
    params:
        inputdir = 'out/demultiplex_qc/{marker}',
        outputdir = 'out/demultiplex_qc/{marker}'
    log:
        'out/demultiplex_qc/{marker}/multiqc_report.log'
    conda:
        'envs/multiqc-1.27.1.yaml'
    shell:
        'multiqc {config[multiqc]} --force -o {params.outputdir} {params.inputdir} 2>{log}'

################################
# remove Ns and check primer orientation in preparation for running cutadapt

rule all_filter_Ns:
    input:
        [os.path.join('out', 'demultiplex_no_Ns', marker, sample + '-' + file) 
            for marker in MARKERS 
            for sample in SAMPLES[marker] 
            for file in ['R1.fastq.gz','R2.fastq.gz','summary.txt']]

rule filter_Ns:
    input:
        read1='out/demultiplex/{marker}/{sample}-R1.fastq.gz',
        read2='out/demultiplex/{marker}/{sample}-R2.fastq.gz'
    output:
        read1='out/demultiplex_no_Ns/{marker}/{sample}-R1.fastq.gz',
        read2='out/demultiplex_no_Ns/{marker}/{sample}-R2.fastq.gz',
        summary='out/demultiplex_no_Ns/{marker}/{sample}-summary.txt'
    params:
        primer_fwd=lambda wildcards: config['primers'][wildcards.marker]['forward'],
        primer_rev=lambda wildcards: config['primers'][wildcards.marker]['reverse']
    conda:
        'envs/dada2-1.34.0.yaml'
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
        [os.path.join('out', 'cutadapt', marker, sample + '-' + file) 
            for marker in MARKERS 
            for sample in SAMPLES[marker] 
            for file in ['R1.fastq.gz','R2.fastq.gz']]

rule cutadapt:
    input:
        read1='out/demultiplex_no_Ns/{marker}/{sample}-R1.fastq.gz',
        read2='out/demultiplex_no_Ns/{marker}/{sample}-R2.fastq.gz'
    output:
        read1='out/cutadapt/{marker}/{sample}-R1.fastq.gz',
        read2='out/cutadapt/{marker}/{sample}-R2.fastq.gz'
    log:
        'out/cutadapt/{marker}/{sample}.log'
    params:
        primer_fwd = lambda wildcards: config['primers'][wildcards.marker]['forward'],
        primer_rev = lambda wildcards: config['primers'][wildcards.marker]['reverse'],
        min_length = config['cutadapt']['min_length']
    conda:
        'envs/cutadapt-5.0.yaml'
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
        [os.path.join('out', 'read_quality_profiles', marker, sample + '.pdf') 
            for marker in MARKERS 
            for sample in SAMPLES[marker]]

rule read_quality_profiles:
    input:
        read1='out/cutadapt/{marker}/{sample}-R1.fastq.gz',
        read2='out/cutadapt/{marker}/{sample}-R2.fastq.gz'
    output:
        'out/read_quality_profiles/{marker}/{sample}.pdf'
    conda:
        'envs/dada2-1.34.0.yaml'
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
        [os.path.join('out', 'filterAndTrim', marker, sample + '-' + read + '.fastq.gz') 
            for marker in MARKERS 
            for sample in SAMPLES[marker] 
            for read in ['R1','R2']]

rule filter_and_trim_16S:
    input:
        read1 = 'out/cutadapt/16S/{sample}-R1.fastq.gz',
        read2 = 'out/cutadapt/16S/{sample}-R2.fastq.gz'
    output:
        read1 = 'out/filterAndTrim/16S/{sample}-R1.fastq.gz',
        read2 = 'out/filterAndTrim/16S/{sample}-R2.fastq.gz'
    params:
        maxN = config['fastqPairedFilter']['16S']['maxN'],
        truncQ = config['fastqPairedFilter']['16S']['truncQ'],
        maxEE_read1 = config['fastqPairedFilter']['16S']['maxEE_read1'],
        maxEE_read2 = config['fastqPairedFilter']['16S']['maxEE_read2'],
        truncLen_read1 = config['fastqPairedFilter']['16S']['truncLen_read1'],
        truncLen_read2 = config['fastqPairedFilter']['16S']['truncLen_read2']
    conda:
        'envs/dada2-1.34.0.yaml'
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
        read1 = 'out/cutadapt/ITS/{sample}-R1.fastq.gz',
        read2 = 'out/cutadapt/ITS/{sample}-R2.fastq.gz'
    output:
        read1 = 'out/filterAndTrim/ITS/{sample}-R1.fastq.gz',
        read2 = 'out/filterAndTrim/ITS/{sample}-R2.fastq.gz'
    params:
        maxN = config['fastqPairedFilter']['ITS']['maxN'],
        truncQ = config['fastqPairedFilter']['ITS']['truncQ'],
        maxEE_read1 = config['fastqPairedFilter']['ITS']['maxEE_read1'],
        maxEE_read2 = config['fastqPairedFilter']['ITS']['maxEE_read2'],
        minLen = config['fastqPairedFilter']['ITS']['minLen']
    conda:
        'envs/dada2-1.34.0.yaml'
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
        [os.path.join('out', 'learnerrors', marker, read, sample + '-' + read + '.fastq.gz') 
            for marker in MARKERS 
            for sample in SAMPLES[marker] 
            for read in ['R1','R2']]

rule learnerrors_inputs:
    input:
        read1 = 'out/filterAndTrim/{marker}/{sample}-R1.fastq.gz',
        read2 = 'out/filterAndTrim/{marker}/{sample}-R2.fastq.gz'
    output:
        read1 = 'out/learnerrors/{marker}/R1/{sample}-R1.fastq.gz',
        read2 = 'out/learnerrors/{marker}/R2/{sample}-R2.fastq.gz'
    shell:
        '''
        ln -s ../../../../{input.read1} {output.read1}
        ln -s ../../../../{input.read2} {output.read2}
        '''

rule all_learnerrors:
    input:
        # one job for each read direction for each marker
        [os.path.join('out', marker, 'learnerrors', read + '_error_profile.rds') 
            for marker in MARKERS 
            for read in ['R1','R2']]

rule learnerrors:
    input:
        # all the fastq files for the specified read direction (R1 or R2) for the specified marker (16S or ITS)
        lambda wildcards: [os.path.join('out', 'learnerrors', wildcards.marker, wildcards.read, sample + '-' + wildcards.read + '.fastq.gz') 
            for sample in SAMPLES[wildcards.marker]]
    output:
        rds = 'out/learnerrors/{marker}/{read}_error_profile.rds',
        pdf = 'out/learnerrors/{marker}/{read}_error_profile.pdf'
    log:
        'out/learnerrors/{marker}/{read}_error_profile.log'
    params:
        input_dir = 'out/learnerrors/{marker}/{read}',
        nbases = config['learnErrors']['nbases']
    conda:
        'envs/dada2-1.34.0.yaml'
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
        [os.path.join('out', 'dada', marker, sample + '.rds') 
            for marker in MARKERS 
            for sample in SAMPLES[marker]]

# this rule appears to be memory intensive, especially for samples with many reads e.g the respiration samples
# -- it may be necessary to limit execution to a couple of jobs at a time to ensure enough memory is available
rule dada:
    input:
        fastq_read1 = 'out/learnerrors/{marker}/R1/{sample}-R1.fastq.gz',
        fastq_read2 = 'out/learnerrors/{marker}/R2/{sample}-R2.fastq.gz',
        error_profile_read1 = 'out/learnerrors/{marker}/R1_error_profile.rds',
        error_profile_read2 = 'out/learnerrors/{marker}/R2_error_profile.rds'
    output:
        'out/dada/{marker}/{sample}.rds'
    log:
        'out/dada/{marker}/{sample}.log'
    conda:
        'envs/dada2-1.34.0.yaml'
    threads:
        config['dada']['threads']
    shell:
        '''
        # test input file and skip processing if empty
        if [ $(zcat {input.fastq_read1} | wc -l) -eq 0 ]
        then
            touch {output}
        else
            Rscript code/dada.R {input.fastq_read1} {input.fastq_read2} \
                {input.error_profile_read1} {input.error_profile_read2} {output} {threads} >{log}
        fi
        '''

rule all_dada_merge:
    input:
        expand('out/dada-merge/{marker}/sequence_table.rds', marker = MARKERS)

rule dada_merge_samples:
    input:
        read1_fastq_files = lambda wildcards: [os.path.join('out', 'learnerrors', wildcards.marker, 'R1', sample + '-R1.fastq.gz') 
            for sample in SAMPLES[wildcards.marker]],
        rds_files = lambda wildcards: [os.path.join('out', 'dada', wildcards.marker, sample + '.rds') 
            for sample in SAMPLES[wildcards.marker]]
    output:
        sample_list = 'out/dada-merge/{marker}/sample_list.txt',
        sequence_table = 'out/dada-merge/{marker}/sequence_table.rds'
    log:
        'out/dada-merge/{marker}/sequence_table.log'
    conda:
        'envs/dada2-1.34.0.yaml'
    shell:
        '''
        # write list of non-empty input files to file
            rm -rf {output.sample_list}
            for READ1_FASTQ_FILE in {input.read1_fastq_files}
            do
                SAMPLE=`basename ${{READ1_FASTQ_FILE%-R1.fastq.gz}}`
                RDS_FILE="out/dada/{wildcards.marker}/${{SAMPLE}}.rds"
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
        expand('out/remove_chimeras/{marker}/sequence_table.rds', marker = MARKERS)

rule remove_chimeras:
    input:
        'out/dada-merge/{marker}/sequence_table.rds'
    output:
        'out/remove_chimeras/{marker}/sequence_table.rds'
    log:
        'out/remove_chimeras/{marker}/sequence_table.log'
    conda:
        'envs/dada2-1.34.0.yaml'
    shell:
        'Rscript code/remove_chimeras.R {input} {output} >{log}'

################################
# assign taxonomy

rule all_taxonomy:
    input:
        expand('out/assign_taxonomy/16S/silva_{matches}.rds', matches = ['single', 'multiple']),
        'out/assign_taxonomy/ITS/unite.rds'

# uses dada2::assignTaxonomy() and dada2::assignSpecies()
# returns a file of unambiguous matches and a file of multiple matches
rule assign_taxonomy_16S:
    input:
        sequence_table = 'out/remove_chimeras/16S/sequence_table.rds',
        assignTaxonomy_refFasta = lambda wildcards: config['databases']['assignTaxonomy'][wildcards.database],
        assignSpecies_refFasta = lambda wildcards: config['databases']['assignSpecies'][wildcards.database]
    output:
        single_match = 'out/assign_taxonomy/16S/{database}_single.rds',
        multiple_matches = 'out/assign_taxonomy/16S/{database}_multiple.rds'
    log:
        'out/assign_taxonomy/16S/{database}.log'
    conda:
        'envs/dada2-1.34.0.yaml'
    threads:
        config['assignTaxonomy']['threads']
    shell:
        '''
        Rscript code/assign_taxonomy_16S.R {input.sequence_table} \
            {input.assignTaxonomy_refFasta} {input.assignSpecies_refFasta} \
            {output.single_match} {output.multiple_matches} {threads} >{log}
        '''

# uses dada2::assignTaxonomy() only
rule assign_taxonomy_ITS:
    input:
        sequence_table = 'out/remove_chimeras/ITS/sequence_table.rds',
        assignTaxonomy_refFasta = lambda wildcards: config['databases']['assignTaxonomy'][wildcards.database]
    output:
        'out/assign_taxonomy/ITS/{database}.rds'
    log:
        'out/assign_taxonomy/ITS/{database}.log'
    conda:
        'envs/dada2-1.34.0.yaml'
    threads:
        config['assignTaxonomy']['threads']
    shell:
        '''
        Rscript code/assign_taxonomy_ITS.R {input.sequence_table} \
            {input.assignTaxonomy_refFasta} {output} {threads} >{log}
        '''

################################
# export as phyloseq object

rule all_export_to_phyloseq:
    input:
        [os.path.join('out', 'phyloseq', marker, 'phyloseq.rds') for marker in MARKERS]

rule export_to_phyloseq_16S:
    input:
        sequence_table = 'out/remove_chimeras/16S/sequence_table.rds',
        metadata = FILEPATHS['16S']['metadata'],
        silva_single = 'out/assign_taxonomy/16S/silva_single.rds',
        silva_multiple = 'out/assign_taxonomy/16S/silva_multiple.rds'
    output:
        phyloseq = 'out/phyloseq/16S/phyloseq.rds',
        fasta = 'out/phyloseq/16S/phyloseq.fasta'
    log:
        'out/phyloseq/16S/phyloseq.log'
    conda:
        'envs/phyloseq-1.50.0.yaml'
    shell:
        '''
        Rscript code/export_to_phyloseq_16S.R {input.sequence_table} {input.metadata} \
            {input.silva_single} {input.silva_multiple} \
            {output.phyloseq} {output.fasta} >{log}
        '''

rule export_to_phyloseq_ITS:
    input:
        sequence_table = 'out/remove_chimeras/ITS/sequence_table.rds',
        metadata = FILEPATHS['ITS']['metadata'],
        unite = 'out/assign_taxonomy/ITS/unite.rds'
    output:
        'out/phyloseq/ITS/phyloseq.rds'
    log:
        'out/phyloseq/ITS/phyloseq.log'
    conda:
        'envs/phyloseq-1.50.0.yaml'
    shell:
        '''
        Rscript code/export_to_phyloseq_ITS.R \
            {input.sequence_table} {input.metadata} {input.unite} \
            {output} >{log}
        '''

################################
rule remove_chloroplasts:
    input:
        rds = 'out/phyloseq/16S/phyloseq.rds'
    output:
        rds = 'out/phyloseq/16S/phyloseq_cleaned.rds',
        fasta = 'out/phyloseq/16S/phyloseq_cleaned.fasta'
    conda:
        'envs/phyloseq-1.50.0.yaml'
    shell:
        'Rscript code/remove_chloroplasts.R {input.rds} {output.rds} {output.fasta}'

################################
# Note: SEPP needs to be installed separately (not through conda) in the `/software/sepp-package` folder.
# The `sepp_env` conda environment provides a suitable runtime environment but does not actually install SEPP.
# The script `scripts/install_sepp.sh` provides suitable commands for installing SEPP.
rule run_sepp:
    input:
        'out/phyloseq/16S/phyloseq_cleaned.fasta'
    output:
        'out/phyloseq/16S/phyloseq_cleaned_placement.tog.tre'
    conda:
        'envs/sepp_env.yaml'
    threads:
        config['sepp']['threads']
    shell:
        # running sepp by cd'ing into this directory seems to be the way to do it but it requires some convoluted relative paths to get it to work here
        '''
        cd out/phyloseq/16S
        ../../../software/sepp-package/run-sepp.sh phyloseq_cleaned.fasta phyloseq_cleaned -x {threads}
        '''

################################
rule combine_data:
    input:
        phyloseq_16S = 'out/phyloseq/16S/phyloseq_cleaned.rds',
        sepp_tree_16S = 'out/phyloseq/16S/phyloseq_cleaned_placement.tog.tre',
        phyloseq_ITS = 'out/phyloseq/ITS/phyloseq.rds'
    output:
        'out/phyloseq/phyloseq_list.rds'
    log:
        'out/phyloseq/phyloseq_list.log'
    conda:
       	'envs/phyloseq-1.50.0.yaml'
    shell:
        'Rscript code/combine_data.R {input.phyloseq_16S} {input.sepp_tree_16S} {input.phyloseq_ITS} {output} >{log}'

################################
# get sequence counts at various stages of the pipeline

# sequence counts prior to demultiplexing

rule sequence_counts_before_demultiplexing:
    input:
        lambda wildcards: FILEPATHS[wildcards.marker]['read1']
    output:
        'out/sequence_counts/before_demultiplexing/{marker}/{marker}.csv'
    shell:
        '''
        STAGE='before demultiplexing'
        MARKER={wildcards.marker}
        SAMPLE=`basename {input}`
        SEQS=`echo $(zcat {input} | wc -l)/4|bc`
        echo ${{STAGE}},${{MARKER}},${{SAMPLE}},${{SEQS}} >{output}
        '''

rule aggregate_sequence_counts_before_demultiplexing:
    input:
        lambda wildcards: [os.path.join('out', 'sequence_counts', 'before_demultiplexing', marker, marker + '.csv') for marker in MARKERS]
    output:
        'out/sequence_counts/before_demultiplexing/aggregated.csv'
    shell:
        '''
        echo 'stage,marker,sample,reads' >{output}
        for f in {input}; do cat $f >>{output}; done
        '''

# sequence counts after demultiplexing

rule sequence_counts_after_demultiplexing:
    input:
        'out/demultiplex/{marker}/{sample}-R1.fastq.gz'
    output:
        'out/sequence_counts/after_demultiplexing/{marker}/{sample}.csv'
    shell:
        '''
        STAGE='after demultiplexing'
        MARKER={wildcards.marker}
        SAMPLE=`basename {input} -R1.fastq.gz`
        SEQS=`echo $(zcat {input} | wc -l)/4|bc`
        echo ${{STAGE}},${{MARKER}},${{SAMPLE}},${{SEQS}} >{output}
        '''

rule aggregate_sequence_counts_after_demultiplexing:
    input:
        lambda wildcards: [os.path.join('out', 'sequence_counts', 'after_demultiplexing', marker, sample + '.csv') 
            for marker in MARKERS 
            for sample in SAMPLES[marker]]
    output:
        'out/sequence_counts/after_demultiplexing/aggregated.csv'
    shell:
        '''
        echo 'stage,marker,sample,reads' >{output}
        for f in {input}; do cat $f >>{output}; done
        '''

# sequence counts after filtering Ns

rule sequence_counts_after_filtering_Ns:
    input:
        'out/demultiplex_no_Ns/{marker}/{sample}-R1.fastq.gz'
    output:
        'out/sequence_counts/after_filtering_Ns/{marker}/{sample}.csv'
    shell:
        '''
        STAGE='after filtering Ns'
        MARKER={wildcards.marker}
        SAMPLE=`basename {input} -R1.fastq.gz`
        SEQS=`echo $(zcat {input} | wc -l)/4|bc`
        echo ${{STAGE}},${{MARKER}},${{SAMPLE}},${{SEQS}} >{output}
        '''

rule aggregate_sequence_counts_after_filtering_Ns:
    input:
        lambda wildcards: [os.path.join('out', 'sequence_counts', 'after_filtering_Ns', marker, sample + '.csv') 
            for marker in MARKERS 
            for sample in SAMPLES[marker]]
    output:
        'out/sequence_counts/after_filtering_Ns/aggregated.csv'
    shell:
        '''
        echo 'marker,stage,sample,reads' >{output}
        for f in {input}; do cat $f >>{output}; done
        '''

# sequence counts after running cutadapt

rule sequence_counts_after_cutadapt:
    input:
        'out/cutadapt/{marker}/{sample}-R1.fastq.gz'
    output:
        'out/sequence_counts/after_cutadapt/{marker}/{sample}.csv'
    shell:
        '''
        STAGE='after cutadapt'
        MARKER={wildcards.marker}
        SAMPLE=`basename {input} -R1.fastq.gz`
        SEQS=`echo $(zcat {input} | wc -l)/4|bc`
        echo ${{STAGE}},${{MARKER}},${{SAMPLE}},${{SEQS}} >{output}
        '''

rule aggregate_sequence_counts_after_cutadapt:
    input:
        lambda wildcards: [os.path.join('out', 'sequence_counts', 'after_cutadapt', marker, sample + '.csv') 
            for marker in MARKERS 
            for sample in SAMPLES[marker]]
    output:
        'out/sequence_counts/after_cutadapt/aggregated.csv'
    shell:
        '''
        echo 'stage,marker,sample,reads' >{output}
        for f in {input}; do cat $f >>{output}; done
        '''

# sequence counts after filterAndTrim

rule sequence_counts_after_filterAndTrim:
    input:
        'out/filterAndTrim/{marker}/{sample}-R1.fastq.gz'
    output:
        'out/sequence_counts/after_filterAndTrim/{marker}/{sample}.csv'
    shell:
        '''
        STAGE='after filterAndTrim'
        MARKER={wildcards.marker}
        SAMPLE=`basename {input} -R1.fastq.gz`
        SEQS=`echo $(zcat {input} | wc -l)/4|bc`
        echo ${{STAGE}},${{MARKER}},${{SAMPLE}},${{SEQS}} >{output}
        '''

rule aggregate_sequence_counts_after_filterAndTrim:
    input:
        lambda wildcards: [os.path.join('out', 'sequence_counts', 'after_filterAndTrim', marker, sample + '.csv') for marker in MARKERS for sample in SAMPLES[marker]]
    output:
        'out/sequence_counts/after_filterAndTrim/aggregated.csv'
    shell:
        '''
        echo 'stage,marker,sample,reads' >{output}
        for f in {input}; do cat $f >>{output}; done
        '''

# sequence counts after running dada

rule sequence_counts_after_dada:
    input:
        'out/dada-merge/{marker}/sequence_table.rds'
    output:
        'out/sequence_counts/after_dada/{marker}/{marker}.csv'
    params:
        stage = 'after dada'
    conda:
        'envs/dada2-1.34.0.yaml'
    shell:
        'Rscript code/sequence_counts_after_dada.R {input} {output} "{params.stage}" {wildcards.marker}'

rule aggregate_sequence_counts_after_dada:
    input:
        expand('out/sequence_counts/after_dada/{marker}/{marker}.csv', marker = MARKERS)
    output:
        'out/sequence_counts/after_dada/aggregated.csv'
    shell:
        '''
        echo 'stage,marker,sample,reads' >{output}
        for f in {input}; do tail -n +2 $f >>{output}; done
        '''

# sequence counts after removing chimeras (same Rscript as after_dada)

rule sequence_counts_after_removing_chimeras:
    input:
        'out/remove_chimeras/{marker}/sequence_table.rds'
    output:
        'out/sequence_counts/after_removing_chimeras/{marker}/{marker}.csv'
    params:
        stage = 'after removing chimeras'
    conda:
        'envs/dada2-1.34.0.yaml'
    shell:
        'Rscript code/sequence_counts_after_dada.R {input} {output} "{params.stage}" {wildcards.marker}'

rule aggregate_sequence_counts_after_removing_chimeras:
    input:
        expand('out/sequence_counts/after_removing_chimeras/{marker}/{marker}.csv', marker = MARKERS)
    output:
        'out/sequence_counts/after_removing_chimeras/aggregated.csv'
    shell:
        '''
        echo 'stage,marker,sample,reads' >{output}
        for f in {input}; do tail -n +2 $f >>{output}; done
        '''

# sequence counts after cleaning

rule sequence_counts_after_cleaning:
    input:
        lambda wildcards: os.path.join('out', 'phyloseq', wildcards.marker, 'phyloseq_cleaned.rds') 
            if wildcards.marker == '16S' 
            else os.path.join('out', 'phyloseq', wildcards.marker, 'phyloseq.rds')
    output:
        'out/sequence_counts/after_cleaning/{marker}/{marker}.csv'
    params:
        stage = 'after cleaning'
    conda:
        'envs/phyloseq-1.50.0.yaml'
    shell:
        'Rscript code/sequence_counts_after_cleaning.R {input} {output} "{params.stage}" {wildcards.marker}'

rule aggregate_sequence_counts_after_cleaning:
    input:
        expand('out/sequence_counts/after_cleaning/{marker}/{marker}.csv', marker = MARKERS)
    output:
        'out/sequence_counts/after_cleaning/aggregated.csv'
    shell:
        '''
        echo 'stage,marker,sample,reads' >{output}
        for f in {input}; do tail -n +2 $f >>{output}; done
        '''

# sequence counts after combining (should be the same as post-cleaning)

rule sequence_counts_after_combining:
    input:
        'out/phyloseq/phyloseq_list.rds'
    output:
        'out/sequence_counts/after_combining/aggregated.csv'
    params:
        stage = 'after combining'
    conda:
        'envs/phyloseq-1.50.0.yaml'
    shell:
        'Rscript code/sequence_counts_after_combining.R {input} {output} "{params.stage}"'

# collate sequence counts

rule aggregate_sequence_counts_all_stages:
    input:
        before_demultiplexing = 'out/sequence_counts/before_demultiplexing/aggregated.csv',
        after_demultiplexing = 'out/sequence_counts/after_demultiplexing/aggregated.csv',
        after_filtering_Ns = 'out/sequence_counts/after_filtering_Ns/aggregated.csv',
        after_cutadapt = 'out/sequence_counts/after_cutadapt/aggregated.csv',
        after_filterAndTrim = 'out/sequence_counts/after_filterAndTrim/aggregated.csv',
        after_dada = 'out/sequence_counts/after_dada/aggregated.csv',
        after_removing_chimeras = 'out/sequence_counts/after_removing_chimeras/aggregated.csv',
        after_cleaning = 'out/sequence_counts/after_cleaning/aggregated.csv',
        after_combining = 'out/sequence_counts/after_combining/aggregated.csv'
    output:
        'out/sequence_counts/aggregated.csv'
    shell:
        '''
        echo 'stage,marker,sample,reads' >{output}
        for f in {input}; do tail -n +2 $f >>{output}; done
        '''

################################
