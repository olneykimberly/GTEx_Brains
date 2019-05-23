# MDS: distances correspond to leading log-fold-changes between each pair of RNA samples.
# Leading log-fold-change is the average (root-mean-square) of the largest absolute log-foldchanges between each pair of samples.

METADATA = file.path('/scratch/mjpete11/GTEx/', 'Metadata.csv')
COUNTS = file.path('/scratch/mjpete11/GTEx/Data_Exploration/', 'Count_Matrix.tsv')
PLOT_DIR = '/scratch/mjpete11/GTEx/Data_Exploration/Sex_Tissue_Age/Multidimensional_Scaling/Plots/'
PLOT_FILE = 'MDS_Tissue_Plots.pdf'

# Load packages                                                                 
library(readr)
library(stringr)
library(dplyr)
library(limma)
library(edgeR)

# Read Metadata CSV.                                                            
samples = read.csv(METADATA, header = TRUE)
samples

# Set rownames of metadata object equal to sample names.                        
rownames(samples) <- samples$Sample                                             

# Read in count matrix
cts <- read.table(COUNTS, sep="\t")

# Replace . to - in colnames
colnames(cts) <- str_replace_all(colnames(cts),pattern = "\\.","-")

# Create DGEList object
y <- DGEList(cts)

# Create every combination of tissue matrices
# Make list of lists of samples for each tissue
tissue_lst <- list()
for(i in 1:length(levels(samples$Tissue))){
  tissue_lst[[i]] <- as.vector(samples[samples$Tissue == levels(samples$Tissue)[i], "Sample"])
}
# Rename lists in list as tissue names
names(tissue_lst) <- levels(samples$Tissue)

# List of metadata for each tissue
meta <- list()
for(i in 1:length(levels(samples$Tissue))){
  meta[[i]] <- samples[samples$Tissue == levels(samples$Tissue)[i],]
}
names(meta) <- levels(samples$Tissue)

# Get log of count matrix and convert to data frame
log_cts <- as.data.frame(cpm(y, log=TRUE))

# Split log count df into seperate dfs classified by tissue
log_count_lst <- list()
for(i in 1:length(tissue_lst)){
  log_count_lst[[i]] <- log_cts[,which(colnames(log_cts) %in% tissue_lst[[i]])]
}

names(log_count_lst) <- names(tissue_lst)

# MDS plots
setwd(PLOT_DIR)
pdf(PLOT_FILE)

# Plot constants 
colors <- c("blue", "darkgreen")
par(mar=c(8, 4.1, 4.1, 8), xpd=TRUE) # margins: bottom, left, top, and right

# MDS plot function for dimensions 1 and 2
MDS_k2_Fun <- function(x, y) {
  plots <- plotMDS(x,
                   top = 1000, 
                   pch = 16, 
                   cex = 1, 
                   dim.plot = c(1,2), 
                   col = colors[samples$Sex], # Assuming colors assigned in alphabetical order; females should be blue
                   main = paste('MDS Plot: Dim 1 and 2; ', y))
           legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors, xpd=TRUE)
}

# MDS plot function for dimensions 3 and 4
MDS_k4_Fun <- function(x, y) {
  plots <- plotMDS(x,
                   top = 1000, 
                   pch = 16, 
                   cex = 1, 
                   dim.plot = c(3,4), 
                   col = colors[samples$Sex], 
                   main = paste('MDS Plot: Dim 3 and 4; ', y))
           legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors, xpd=TRUE)
}

# Apply plot function to list of log count dfs
Map(MDS_k2_Fun, x = log_count_lst, y = names(log_count_lst))

Map(MDS_k4_Fun, x = log_count_lst, y = names(log_count_lst))

dev.off()
