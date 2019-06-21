# Create count matrix from hisat + stringtie results. 
# tximport will compute counts from the coverage information, by reversing the formula that StringTie uses to calculate coverage.

library(tximport)                                                               
library(GenomicFeatures) 
library(AnnotationDbi)
library(edgeR)   
library(readr) 

# Read Metadata CSV.                                                            
samples = read.csv(file.path("/scratch/mjpete11/GTEx/Metadata/", "Metadata.csv"), header = TRUE)

# Drop sample that doesn't exist in config, but does in metadata. 
samples <- samples[-c(53),]  # 53: GTEX-ZAB4-0011-R4a-SM-4SOKB: Amygdala, male

# Females and males were aligned to different refernce genomes
# Try to load as seperate objects: https://support.bioconductor.org/p/96588/

fem_samples <- samples[which(samples$Sex == "Female"),]
male_samples <- samples[which(samples$Sex == "Male"),]

# Set rownames of metadata object equal to sample names.                        
rownames(fem_samples) <- fem_samples$Sample 
rownames(male_samples) <- male_samples$Sample

# Set path to quant.sf files and check that all files are there.                
fem_files = file.path("/scratch/mjpete11/GTEx/Stringtie_Quants", fem_samples$Sample, "t_data.ctab")
all(file.exists(fem_files)) 

male_files = file.path("/scratch/mjpete11/GTEx/Stringtie_Quants", male_samples$Sample, "t_data.ctab")
all(file.exists(male_files)) 

# The tx2gene table should connect transcripts to genes, and can be pulled out of one of the t_data.ctab files
fem_tmp <- read_tsv(fem_files[1])
male_tmp <- read_tsv(male_files[1])

fem_tx2gene <- fem_tmp[, c("t_name", "gene_name")]
male_tx2gene <- male_tmp[, c("t_name", "gene_name")]

fem_txi <- tximport(fem_files, type = "stringtie", tx2gene = tx2gene)
male_txi <- tximport(male_files, type = "stringtie", tx2gene = tx2gene)

# Get counts
fem_cts <- fem_txi$counts
male_cts <- male_txi$counts

df_fem_cts <- data.frame(fem_cts, check.names = FALSE)
df_male_cts <- data.frame(male_cts, check.names = FALSE)

# Combine to one df
df_cts <- Merge(df_fem_cts, df_male_cts, by=NULL)

# Write to TSV
write.table(df_cts, file = "Stringtie_Count_Matrix.tsv", sep = "\t", row.names = TRUE, quote = FALSE)
