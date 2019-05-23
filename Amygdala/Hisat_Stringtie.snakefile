configfile: "/scratch/mjpete11/GTEx/Configs/GTEx_Brain_v1.8.config"

import os

# TOOLS
# Variables should point to full paths to tools.
# Currently assumes all are available in one's PATH
# (e.g., installed with Bioconda)
HISAT2_BUILD_PATH = "hisat2-build"
HISAT2_PATH = "hisat2"
SAMTOOLS_PATH = "samtools"
XYALIGN_ENV_PATH = "xyalign_env"
XYALIGN_PATH = "xyalign"
STRINGTIE_PATH = "stringtie"

# Directory variables
FASTQ_DIRECTORY = "/mnt/storage/public/dbgap-8834/brain_amygdala/"

# Samples
XX_SAMPLES = config["Brain_Amygdala_Female_RNA"]
XY_SAMPLES_wo = config["Brain_Amygdala_Male_RNA_wo_ZAB4"]
XY_SAMPLES = config["Brain_Amygdala_Male_RNA"]
SAMPLES = XX_SAMPLES + XY_SAMPLES_wo


rule all:
  input:
    "multiqc_results/multiqc_report.html",
    "multiqc_trimmed_results/multiqc_report.html",
        expand("stringtie_results/sample_{sample}/{sample}_{assembly}.assembled_transcripts.secondpass.gtf", assembly=["hg38"], sample=SAMPLES)

rule xyalign_prepare_ref:
    input:
        lambda wildcards: config["ref_genome"][wildcards.assembly]
    output:
        xx = "xyalign/reference/{assembly}_xx.fa",
        xy = "xyalign/reference/{assembly}_xy.fa"
    params:
        xyalign_env = XYALIGN_ENV_PATH,
        xyalign = XYALIGN_PATH
    shell:
        "source activate {params.xyalign_env} && "
        "xyalign --PREPARE_REFERENCE --ref {input} --output_dir xyalign "
        "--reference_mask hg38_PAR_mask.bed --x_chromosome chrX --y_chromosome chrY "
        "--xx_ref_out {wildcards.assembly}_xx.fa "
        "--xy_ref_out {wildcards.assembly}_xy.fa"

rule hisat2_index:
    input:
        xx = "xyalign/reference/{assembly}_xx.fa",
        xy = "xyalign/reference/{assembly}_xy.fa"
    output:
        xx = "hisat2_index/{assembly}_xx.6.ht2",
        xy = "hisat2_index/{assembly}_xy.6.ht2"
    params:
        hisat2 = HISAT2_BUILD_PATH,
        base_name1 = "hisat2_index/{assembly}_xx",
        base_name2 = "hisat2_index/{assembly}_xy"
    run:
        shell("{params.hisat2} {input.xx} {params.base_name1}")
        shell("{params.hisat2} {input.xy} {params.base_name2}")

rule hisat2_align_reads:
    input:
        paired_1 = "trimmed_fastqs/{sample}_trimmomatic_trimmed_paired_1.fastq.gz",
        paired_2 = "trimmed_fastqs/{sample}_trimmomatic_trimmed_paired_2.fastq.gz",
        xx = "hisat2_index/{assembly}_xx.6.ht2",
        xy = "hisat2_index/{assembly}_xy.6.ht2"
    output:
        "bams/{sample}_{assembly}.sorted.bam"
    params:
        hisat2 = HISAT2_PATH,
        threads = 4,
        base_name_xx = "hisat2_index/{assembly}_xx",
        base_name_xy = "hisat2_index/{assembly}_xy",
        samtools = SAMTOOLS_PATH
        Threads: 4
    run:
        if wildcards.sample in XY_SAMPLES_wo:
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

rule stringtie_first_pass:
    input:
        bam = "bams/{sample}_{assembly}.sorted.bam",
        gff = "reference/{assembly}.gff"
    output:
        "stringtie_results/sample_{sample}/{sample}_{assembly}.assembled_transcripts.firstpass.gtf"
    threads: 4 # what is this?
    params:
        stringtie = STRINGTIE_PATH,
        threads = 4
    shell:
        "{params.stringtie} {input.bam} -o {output} -p {params.threads} "
        "-G {input.gff}"

rule create_stringtie_merged_list:
    input:
        lambda wildcards: expand(
            "stringtie_results/sample_{sample}/{sample}_{assembly}.assembled_transcripts.firstpass.gtf",
            assembly = WILDCARDS.GENOME,
            sample = SAMPLES)
    output:
        "stringtie_results/{genome}_gtflist.txt"
    run:
        shell("echo -n > {output}")
        for i in input:
            shell("echo {} >> {{output}}".format(i))

rule stringtie_merge:
    input:
        stringtie_list = "stringtie_results/{genome}_gtflist.txt",
        gff = "reference/{genome}.gff"
    output:
        "stringtie_results/{genome}.merged.gtf"
    threads: 4
    params:
        stringtie = STRINGTIE_PATH,
        threads = 4
    shell:
        "{params.stringtie} --merge {input.stringtie_list} -o {output} "
        "-p {params.threads} -G {input.gff}"

rule stringtie_second_pass:
    input:
        bam = "bams/{sample}_{assembly}.sorted.bam",
        gff = "stringtie_results/{assembly}.merged.gtf"
    output:
        assembled_transcripts = "stringtie_results/sample_{sample}/{sample}_{assembly}.assembled_transcripts.secondpass.gtf",
        gene_abundances = "stringtie_results/sample_{sample}/{sample}_{assembly}.gene_abundances.secondpass.txt",
        fully_covered_transcripts = "stringtie_results/sample_{sample}/{sample}_{assembly}.fully_covered_transcripts.secondpass.gtf"
    threads: 4
    params:
        stringtie = STRINGTIE_PATH,
        threads = 4
    shell:
        "{params.stringtie} {input.bam} -p {params.threads} "
        "-G {input.gff} -B -e "
        "-o {output.assembled_transcripts} "
        "-A {output.gene_abundances} "
        "-C {output.fully_covered_transcripts}"
