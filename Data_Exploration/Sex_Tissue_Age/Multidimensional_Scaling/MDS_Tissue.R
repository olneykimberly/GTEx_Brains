# This script is for plotting MDS dimensions 1:4.

METADATA = file.path('/scratch/mjpete11/GTEx/', 'Metadata.csv')
COUNTS = file.path('/scratch/mjpete11/GTEx/Data_Exploration/', 'Count_Matrix.tsv')
PLOT_DIR = 'Plots/'
PLOT_FILE = 'Tissue_MDS_Sex_RD.pdf'

# Load packages                                                                 
library(readr)
library(stringr)
library(dplyr)
library(limma)
library(edgeR)

# Read Metadata CSV.                                                            
samples = read.csv(METADATA, header = TRUE)
samples                                                                         

# Create every combination of tissue matrices
# Make list of lists of samples for each tissue
tissue_lst <- list()
for(i in 1:length(levels(samples$Tissue))){
  tissue_lst[[i]] <- as.vector(samples[samples$Tissue == levels(samples$Tissue)[i], "Sample"])
}
# Rename lists in list as tissue names
names(tissue_lst) <- levels(samples$Tissue)

# Read in counts 
cts <- read.csv(COUNTS, sep = "\t")

# Replace . to - in colnames
colnames(cts) <- str_replace_all(colnames(cts),pattern = "\\.","-")

# List of metadata for each tissue
meta <- list()
for(i in 1:length(levels(samples$Tissue))){
  meta[[i]] <- samples[samples$Tissue == levels(samples$Tissue)[i],]
}
names(meta) <- levels(samples$Tissue)

# Make list of data frames of counts from samples per tissue
count_lst <- list()
for(i in 1:length(tissue_lst)){
  count_lst[[i]] <- cts[,which(colnames(cts) %in% tissue_lst[[i]])]
}

names(count_lst) <- names(tissue_lst)

# MDS plots
setwd(PLOT_DIR)
pdf(PLOT_FILE)

# Plot

# Amygdala dim 1 & 2
colors <- c("blue", "darkgreen")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) # Add extra space to right of plot area; change clipping to figure
plotMDS(count_lst["Amygdala"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(1,2), 
        col = colors[samples$Sex], # Assuming colors assigned in alphabetical order; females should be blue
        main = 'MDS plot dim 1 and 2: Amygdala')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Amygdala dim 3 & 4
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) # Add extra space to right of plot area; change clipping to figure
plotMDS(count_lst["Amygdala"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(3,4), 
        col = colors[samples$Sex], # Assuming colors assigned in alphabetical order; females should be blue
        main = 'MDS plot dim 3 and 4: Amygdala')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Anterior dim 1 & 2
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotMDS(count_lst["Anterior"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(1,2), 
        col = colors[samples$Sex],
        main = 'MDS plot dim 1 and 2: Anterior')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Anterior dim 3 & 4
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) 
plotMDS(count_lst["Anterior"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(3,4), 
        col = colors[samples$Sex],
        main = 'MDS plot dim 3 and 4: Anterior')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Caudate dim 1 & 2
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotMDS(count_lst["Caudate"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(1,2), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 1 and 2: Caudate')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Caudate dim 3 & 4
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) 
plotMDS(count_lst["Caudate"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(3,4), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 3 and 4: Caudate')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Frontal Cortex dim 1 & 2
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotMDS(count_lst["Frontal_Cortex"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(1,2), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 1 and 2: Frontal Cortex')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Frontal Cortex dim 3 & 4
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) 
plotMDS(count_lst["Frontal_Cortex"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(3,4), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 3 and 4: Frontal Cortex')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Cortex dim 1 & 2
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotMDS(count_lst["Cortex"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(1,2), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 1 and 2: Cortex')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Cortex dim 3 & 4
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) 
plotMDS(count_lst["Cortex"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(3,4), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 3 and 4: Cortex')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Cerebellar dim 1 & 2
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotMDS(count_lst["Cerebellar"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(1,2), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 1 and 2: Cerebellar')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Cerebellar dim 3 & 4
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) 
plotMDS(count_lst["Cerebellar"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(3,4), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 3 and 4: Cerebellar')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Cerebellum dim 1 & 2
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotMDS(count_lst["Cerebellum"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(1,2), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 1 and 2: Cerebellum')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Cerebellum dim 3 & 4
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) 
plotMDS(count_lst["Cerebellum"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(3,4), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 3 and 4: Cerebellum')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Hypothalamus dim 1 & 2
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotMDS(count_lst["Hypothalamus"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(1,2), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 1 and 2: Hypothalamus')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Hypothalamus dim 3 & 4
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) 
plotMDS(count_lst["Hypothalamus"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(3,4), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 3 and 4: Hypothalamus')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Hippocampus dim 1 & 2
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotMDS(count_lst["Hippocampus"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(1,2), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 1 and 2: Hippocampus')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Hippocampus dim 3 & 4
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) 
plotMDS(count_lst["Hippocampus"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(3,4), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 3 and 4: Hippocampus')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Substantia Nigra dim 1 & 2
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotMDS(count_lst["Substantia_Nigra"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(1,2), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 1 and 2: Substantia Nigra')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Substantia Nigra dim 3 & 4
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) 
plotMDS(count_lst["Substantia_Nigra"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(3,4), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 3 and 4: Substantia Nigra')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Spinal Cord dim 1 & 2
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotMDS(count_lst["Spinal_Cord"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(1,2), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 1 and 2: Spinal Cord')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Spinal Cord dim 3 & 4
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) 
plotMDS(count_lst["Spinal_Cord"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(3,4), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 3 and 4: Spinal Cord')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Nucleus Accumbens dim 1 & 2
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotMDS(count_lst["Nucleus_Accumbens"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(1,2), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 1 and 2: Nucleus Accumbens')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Nucleus Accumbens dim 3 & 4
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) 
plotMDS(count_lst["Nucleus_Accumbens"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(3,4), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 3 and 4: Nucleus Accumbens')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Putamen dim 1 & 2
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotMDS(count_lst["Putamen"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(1,2), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 1 and 2: Putamen')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

# Putamen dim 3 & 4
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) 
plotMDS(count_lst["Putamen"][[1]], 
        top = 1000, 
        pch = 16, 
        cex = 1, 
        dim.plot = c(3,4), 
        col = colors[samples$Sex], 
        main = 'MDS plot dim 3 and 4: Putamen')
legend("topright", inset=c(-0.2,0), legend=levels(samples$Sex), pch=16, cex=1,col=colors)

dev.off()
