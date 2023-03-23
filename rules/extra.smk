rule fastqc:
    input:
        'data/{sample}_{read}.fastq.gz'
    output:
        'out/fastqc/{sample}_{read}_fastqc.html',
        'out/fastqc/{sample}_{read}_fastqc.zip'
    conda:
        'envs/fastqc-0.11.9.yaml'
    shell:
        'fastqc -o out/fastqc {input}'
