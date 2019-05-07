config = ["Male_Spinal_Cord.config.json"]
# Contains dictionary mapping GTEX IDs to SRA IDs

import os

# Directories
INPUT_DIRECTORY = "/data/storage/public/dbgap-8834/malebrain/"
OUTPUT_DIRECTORY = "/data/storage/public/dbgap-8834/brain_spinal_cord/MALE" 

# List of male spinal cord GTEX IDs. 
SPINAL_CORD_MALE = config["Brain_Spinal_Cord_RNA_Male"]

rule all:
    input: 
        expand(
            "{GTEX_ID}/{GTEX_ID}_1.fastq", 
            "{GTEX_ID}/{GTEX_ID}_2.fastq", GTEX_ID=SPINAL_CORD_MALE)


rule make_symlinks:
    input:
        fq1 = lambda wildcards: os.path.join(INPUT_DIRECTORY, [wildcards.sample]["fastq1"]),
        fq2 = lambda wildcards: os.path.join(INPUT_DIRECTORY, [wildcards.sample]["fastq2"])
    output:
        ofq1 = os.path.join(OUTPUT_DIRECTORY, "{GTEX_ID}/{GTEX_ID}_1.fastq"),
        ofq2 = os.path.join(OUTPUT_DIRECTORY, "{GTEX_ID}/{GTEX_ID}_2.fastq")
    shell:
        """
        ln -s {input.fq1} {output.ofq1} && touch -h {output.ofq1};
        ln -s {input.fq2} {output.ofq2} && touch -h {output.ofq2};
        """

