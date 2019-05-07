# This script is to create the count matrix of all the samples.

# Load packages                                                                 
library(tximport)                                                               
library(GenomicFeatures) 
library(AnnotationDbi)
library(edgeR)   
library(readr)

# Read Metadata CSV.                                                            
samples = read.csv(file.path("/scratch/mjpete11/GTEx/", "No_Dup_Metadata.csv"), header = TRUE)

# Set rownames of metadata object equal to sample names.                        
rownames(samples) <- samples$Sample                                             

# Set path to quant.sf files and check that all files are there.                
files = file.path("/scratch/mjpete11/GTEx/All_Salmon_Quants", samples$Sample, "quant.sf")
head(files)                                                                     
all(file.exists(files)) 

# GTEX-13N2G-0011-R2a-SM-5MR4Q (male: substantia nigra) is missing the quant file. 
# Removed from No_Dup_Metadata.csv

# Set the names of the file paths object equal to the sample names.             
names(files) <- samples$Sample                                                  

# Make tx2gene table to map gene IDs to transcript IDs.                         
TxDb <- makeTxDbFromGFF(file = "/scratch/mjpete11/GTEx/gencode.v29.annotation.gtf")
k <- keys(TxDb, keytype = "TXNAME")                                             
tx2gene <- select(TxDb, k, "GENEID", "TXNAME")                                  
head(tx2gene)                                                                   

# Import salmon output files. Set txOut=True to keep transcript expression counts, not gene counts.
# If you want the list of genes, add: txOut = FALSE. To get list of isoforms, add: txOut = TRUE
Txi <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = FALSE)      

cts <- Txi$counts

df_cts <- data.frame(cts, check.names = FALSE)

# Write to TSV
write.table(df_cts, file = "Count_Matrix.tsv", sep = "\t", row.names = TRUE, quote = FALSE)
