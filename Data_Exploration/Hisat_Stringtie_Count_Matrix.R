# Create count matrix from hisat + stringtie results. 
# tximport will compute counts from the coverage information, by reversing the formula that StringTie uses to calculate coverage.

library(tximport)                                                               
library(GenomicFeatures) 
library(AnnotationDbi)
library(edgeR)   
library(readr) 

# Read Metadata CSV.                                                            
samples = read.csv(file.path("/scratch/mjpete11/GTEx/Metadata/", "Metadata.csv"), header = TRUE)

# Set rownames of metadata object equal to sample names.                        
rownames(samples) <- samples$Sample                                             

# Set path to quant.sf files and check that all files are there.                
files = file.path("/scratch/mjpete11/GTEx/Anterior/Hisat_Stringtie/stringtie_results", samples$Sample, "t_data.ctab")
head(files)                                                                     
all(file.exists(files)) 

# The tx2gene table should connect transcripts to genes, and can be pulled out of one of the t_data.ctab files
tmp <- read_tsv(files[1])
tx2gene <- tmp[, c("t_name", "gene_name")]
txi <- tximport(files, type = "stringtie", tx2gene = tx2gene)

# Import stringtie output files.
Txi <- tximport(files, type = "stringtie", tx2gene = tx2gene, txOut = FALSE)      

cts <- Txi$counts

df_cts <- data.frame(cts, check.names = FALSE)

# Write to TSV
write.table(df_cts, file = "Stringtie_Count_Matrix.tsv", sep = "\t", row.names = TRUE, quote = FALSE)
