configfile: "/scratch/mjpete11/GTEx/Configs/GTEx_Brain_v1.8.config"

import os

# TOOLS
# Variables should point to full paths to tools.
# Currently assumes all are available in one's PATH
# (e.g., installed with Bioconda)
HISAT2_PATH = "hisat2"
SAMTOOLS_PATH = "samtools"
STRINGTIE_PATH = "stringtie"

# Directory variables
FASTQ_DIRECTORY = "/mnt/storage/public/dbgap-8834/brain_amygdala/"

# Samples
XX_SAMPLES = config["Brain_Amygdala_RNA_Female"]
XY_SAMPLES = config["Brain_Amygdala_RNA_Male"]
SAMPLES = XX_SAMPLES + XY_SAMPLES

rule all:
    input:
        "bams/{sample}_{assembly}.sorted.bam" 

rule hisat2_align_reads:
    input:
        paired_1 = "trimmed_fastqs/{sample}_trimmomatic_trimmed_paired_1.fastq.gz",
        paired_2 = "trimmed_fastqs/{sample}_trimmomatic_trimmed_paired_2.fastq.gz",
        xx = "/mnt/storage/SAYRES/REFERENCE_GENOMES/GENCODE/GRCh38.p12.genome.XXonly/GRCh38.p12.genome.XXonly.8.ht2",
        xy = "/mnt/storage/SAYRES/REFERENCE_GENOMES/GENCODE/GRCh38.p12.genome.XY/GRCh38.p12.genome.XY.8.ht2"
    output:
        "bams/{sample}_{assembly}.sorted.bam"
    params:
        hisat2 = HISAT2_PATH,
        samtools = SAMTOOLS_PATH,
        threads = 4,
        base_name_xx = "Ch38.p12.genome.XXonly",
        base_name_xy = "GRCh38.p12.genome.XY",
    run:
        if wildcards.sample in XY_SAMPLES:
            shell(
                "{params.hisat2} -p {params.threads} -x {params.base_name_xy} "
                "-1 {input.paired_1} -2 {input.paired_2} | "
                "{params.samtools} view -b - | "
                "{params.samtools} sort -O bam -o {output} -")
        else:
            shell(
                "{params.hisat2} -p {params.threads} -x {params.base_name_xx} "
                "-1 {input.paired_1} -2 {input.paired_2} | "
                "{params.samtools} view -b - | "
                "{params.samtools} sort -O bam -o {output} -")

