# This script looks at differential gene expression between each brain tissue type.

METADATA <- "/scratch/mjpete11/GTEx/Metadata.csv"
COUNT_MATRIX <- "/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrix.tsv"
PLOT_DIR =  "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Tissue/Gene/Unnormalized/Exact_Test/Plots/"

# Load packages                                                                 
library(tximport)                                                               
library(GenomicFeatures)                                                        
library(edgeR) 
library(readr)
library(stringr)
library(gridExtra)
library(grid)

# Read Metadata CSV.                                                            
samples = read.csv(METADATA, header = TRUE)

# Set rownames of metadata object equal to sample names.                        
rownames(samples) <- samples$Sample      

# Read in counts 
cts <- read.csv(COUNT_MATRIX, sep = "\t")

# Replace . to - in colnames
colnames(cts) <- str_replace_all(colnames(cts),pattern = "\\.","-")

# Create design matrix: Step 1
# 26 combos of tissue and sex
Group <- factor(paste(samples$Tissue, samples$Sex, sep="."))
cbind(samples,group=Group)

# Create DGEList object: Tissue only
y <- DGEList(cts, group=samples$Tissue)

#Create deisgn matrix: Step 2
design <- model.matrix(~0+group, data=y$samples) # No intercept
colnames(design) <- levels(y$samples$group)
design

# Filter out lowly expressed genes. (< 7 counts)
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]

# Estimate dispersion
y <- estimateDisp(y, design, robust=TRUE)
sqrt(y$common.dispersion)

# Make contrasts: Tissue vs Tissue
my.contrasts <- makeContrasts(Am.Vs.At = Amygdala - Anterior, 
                              Am.Vs.Ca = Amygdala - Caudate,
                              Am.Vs.Ce = Amygdala - Cerebellar,
                              Am.Vs.Cm = Amygdala - Cerebellum,
                              Am.Vs.Co = Amygdala - Cortex,
                              Am.Vs.Fc = Amygdala - Frontal_Cortex,
                              Am.Vs.Hp = Amygdala - Hippocampus,
                              Am.Vs.Hy = Amygdala - Hypothalamus,
                              Am.Vs.Nc = Amygdala - Nucleus_Accumbens,
                              Am.Vs.Pu = Amygdala - Putamen,
                              Am.Vs.Sp = Amygdala - Spinal_Cord,
                              Am.Vs.Sn = Amygdala - Substantia_Nigra,
                              At.Vs.Ca = Anterior - Caudate,
                              At.Vs.Ce = Anterior - Cerebellar,
                              At.Vs.Cm = Anterior - Cerebellum,
                              At.Vs.Co = Anterior - Cortex,
                              At.Vs.Fc = Anterior - Frontal_Cortex,
                              At.Vs.Hp = Anterior - Hippocampus,
                              At.Vs.Hy = Anterior - Hypothalamus,
                              At.Vs.Nc = Anterior - Nucleus_Accumbens,
                              At.Vs.Pu = Anterior - Putamen,
                              At.Vs.Sp = Anterior - Spinal_Cord,
                              At.Vs.Sn = Anterior - Substantia_Nigra,
                              Ca.Vs.Ce = Caudate - Cerebellar,
                              Ca.Vs.Cm = Caudate - Cerebellum,
                              Ca.Vs.Co = Caudate - Cortex,
                              Ca.Vs.Fc = Caudate - Frontal_Cortex,
                              Ca.Vs.Hp = Caudate - Hippocampus,
                              Ca.Vs.Hy = Caudate - Hypothalamus,
                              Ca.Vs.Na = Caudate - Nucleus_Accumbens,
                              Ca.Vs.Pu = Caudate - Putamen,
                              Ca.Vs.Sp = Caudate - Spinal_Cord,
                              Ca.Vs.Sn = Caudate - Substantia_Nigra,
                              Ce.Vs.Cm = Cerebellar - Cerebellum,
                              Ce.Vs.Co = Cerebellar - Cortex,
                              Ce.Vs.Fc = Cerebellar - Frontal_Cortex,
                              Ce.Vs.Hp = Cerebellar - Hippocampus,
                              Ce.Vs.Hy = Cerebellar - Hypothalamus,
                              Ce.Vs.Na = Cerebellar - Nucleus_Accumbens,
                              Ce.Vs.Pu = Cerebellar - Putamen,
                              Ce.Vs.Sp = Cerebellar - Spinal_Cord,
                              Ce.Vs.Sn = Cerebellar - Substantia_Nigra,
                              Co.Vs.Fc = Cortex - Frontal_Cortex,
                              Co.Vs.Hp = Cortex - Hippocampus,
                              Co.Vs.Hy = Cortex - Hypothalamus,
                              Co.Vs.Na = Cortex - Nucleus_Accumbens,
                              Co.Vs.Pu = Cortex - Putamen,
                              Co.Vs.Sp = Cortex - Spinal_Cord,
                              Co.Vs.Sn = Cortex - Substantia_Nigra,
                              Fc.Vs.Hp = Frontal_Cortex - Hippocampus,
                              Fc.Vs.Hy = Frontal_Cortex - Hypothalamus,
                              Fc.Vs.Na = Frontal_Cortex - Nucleus_Accumbens,
                              Fc.Vs.Pu = Frontal_Cortex - Putamen,
                              Fc.Vs.Sp = Frontal_Cortex - Spinal_Cord,
                              Fc.Vs.Sn = Frontal_Cortex - Substantia_Nigra,
                              Hp.Vs.Hy = Hippocampus - Hypothalamus,
                              Hp.Vs.Na = Hippocampus - Nucleus_Accumbens,
                              Hp.Vs.Pu = Hippocampus - Putamen,
                              Hp.Vs.Sp = Hippocampus - Spinal_Cord,
                              Hp.Vs.Sn = Hippocampus - Substantia_Nigra,
                              Hy.Vs.Na = Hypothalamus - Nucleus_Accumbens,
                              Hy.Vs.Pu = Hypothalamus - Putamen,
                              Hy.Vs.Sp = Hypothalamus - Spinal_Cord,
                              Hy.Vs.Sn = Hypothalamus - Substantia_Nigra,
                              Na.Vs.Pu = Nucleus_Accumbens - Putamen,
                              Na.Vs.Sp = Nucleus_Accumbens - Spinal_Cord,
                              Na.Vs.Sn = Nucleus_Accumbens - Substantia_Nigra,
                              Pu.Vs.Sp = Putamen - Spinal_Cord,
                              Pu.Vs.Sn = Putamen - Substantia_Nigra,
                              Sp.Vs.Sn = Spinal_Cord - Substantia_Nigra,
                              levels=design)

# Test for DGX with Exact Test

# Plot
setwd(PLOT_DIR)

pdf('Unnormalized_Exact_Test.pdf')

# Test for DGX using quasi-liklihood method.
# Amygdala vs Anterior
et.Am.Vs.At <- exactTest(y, c("Anterior", "Amygdala"))
df_Am.Vs.At <- summary(decideTests(et.Am.Vs.At))
grid.newpage()
grid.table(df_Am.Vs.At)
plotMD(et.Am.Vs.At)

# Volcano plot
volcanoData <- cbind(et.Am.Vs.At$table$logFC, -log10(et.Am.Vs.At$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Anterior')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Am.Vs.At$table[,"PValue"], breaks=50, main="Amygdala vs Anterior p-value frequency histogram")

# Amygdala vs Caduate
et.Am.Vs.Ca <- exactTest(y, c("Caudate", "Amygdala"))
df_Am.Vs.Ca <- summary(decideTests(et.Am.Vs.Ca))
grid.newpage()
grid.table(df_Am.Vs.Ca)
plotMD(et.Am.Vs.Ca)

# Volcano plot
volcanoData <- cbind(et.Am.Vs.Ca$table$logFC, -log10(et.Am.Vs.Ca$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Caudate')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Am.Vs.Ca$table[,"PValue"], breaks=50, main="Amygdala vs Caudate p-value frequency histogram")

# Amygdala vs Cerebellar
et.Am.Vs.Ce <- exactTest(y, c("Cerebellar", "Amygdala"))
df_Am.Vs.Ce <- summary(decideTests(et.Am.Vs.Ce))
grid.newpage()
grid.table(df_Am.Vs.Ce)
plotMD(et.Am.Vs.Ce)

# Volcano plot
volcanoData <- cbind(et.Am.Vs.Ce$table$logFC, -log10(et.Am.Vs.Ce$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Cerebellar')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Am.Vs.Ce$table[,"PValue"], breaks=50, main="Amygdala vs Cerebellar p-value frequency histogram")

# Amygdala vs Cerebellum
et.Am.Vs.Cm <- exactTest(y, c("Cerebellum", "Amygdala"))
df_Am.Vs.Cm <- summary(decideTests(et.Am.Vs.Cm))
grid.newpage()
grid.table(df_Am.Vs.Cm)
plotMD(et.Am.Vs.Cm)

# Volcano plot
volcanoData <- cbind(et.Am.Vs.Cm$table$logFC, -log10(et.Am.Vs.Cm$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Cerebellum')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Am.Vs.Cm$table[,"PValue"], breaks=50, main="Amygdala vs Cerebellum p-value frequency histogram")

# Amygdala vs Cortex
et.Am.Vs.Co <- exactTest(y, c("Cortex", "Amygdala"))
df_Am.Vs.Co <- summary(decideTests(et.Am.Vs.Co))
grid.newpage()
grid.table(df_Am.Vs.Co)
plotMD(et.Am.Vs.Co)

# Volcano plot
volcanoData <- cbind(et.Am.Vs.Co$table$logFC, -log10(et.Am.Vs.Co$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Am.Vs.Co$table[,"PValue"], breaks=50, main="Amygdala vs Cortex p-value frequency histogram")

# Amygdala vs Frontal Cortex
et.Am.Vs.Fc <- exactTest(y, c("Frontal_Cortex", "Amygdala"))
df_Am.Vs.Fc <- summary(decideTests(et.Am.Vs.Fc))
grid.newpage()
grid.table(df_Am.Vs.Fc)
plotMD(et.Am.Vs.Fc)

# Volcano plot
volcanoData <- cbind(et.Am.Vs.Fc$table$logFC, -log10(et.Am.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Frontal_Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Am.Vs.Fc$table[,"PValue"], breaks=50, main="Amygdala vs Frontal Cortex p-value frequency histogram")

# Amygdala vs Hippocampus
et.Am.Vs.Hp <- exactTest(y, c("Hippocampus", "Amygdala"))
summary(decideTests(et.Am.Vs.Hp))
plotMD(et.Am.Vs.Hp)

# Volcano plot
volcanoData <- cbind(et.Am.Vs.Hp$table$logFC, -log10(et.Am.Vs.Hp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Am.Vs.Hp$table[,"PValue"], breaks=50, main="Amygdala vs Hippocampus p-value frequency histogram")

# Amygdala vs Hypothalamus
et.Am.Vs.Hy <- exactTest(y, c("Hypothalamus", "Amygdala"))
df_Am.Vs.Hy <- summary(decideTests(et.Am.Vs.Hy))
grid.newpage()
grid.table(df_Am.Vs.Hy)
plotMD(et.Am.Vs.Hy)

# Volcano plot
volcanoData <- cbind(et.Am.Vs.Hy$table$logFC, -log10(et.Am.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Am.Vs.Hy$table[,"PValue"], breaks=50, main="Amygdala vs Hypothalamus p-value frequency histogram")

# Amygdala vs Nucleus Accumbens
et.Am.Vs.Nc <- exactTest(y, c("Nucleus_Accumbens", "Amygdala"))
df_Am.Vs.Nc <- summary(decideTests(et.Am.Vs.Nc))
grid.newpage()
grid.table(df_Am.Vs.Nc)
plotMD(et.Am.Vs.Nc)

# Volcano plot
volcanoData <- cbind(et.Am.Vs.Nc$table$logFC, -log10(et.Am.Vs.Nc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Am.Vs.Nc$table[,"PValue"], breaks=50, main="Amygdala vs Nucleus Accumbens p-value frequency histogram")

# Amygdala vs Putamen
et.Am.Vs.Pu <- exactTest(y, c("Putamen", "Amygdala"))
df_Am.Vs.Pu <- summary(decideTests(et.Am.Vs.Pu))
grid.newpage()
grid.table(df_Am.Vs.Pu)
plotMD(et.Am.Vs.Pu)

# Volcano plot
volcanoData <- cbind(et.Am.Vs.Pu$table$logFC, -log10(et.Am.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Am.Vs.Pu$table[,"PValue"], breaks=50, main="Amygdala vs Putamen p-value frequency histogram")

# Amygdala vs Spinal Cord
et.Am.Vs.Sp <- exactTest(y, c("Spinal_Cord", "Amygdala"))
df_Am.Vs.Sp <- summary(decideTests(et.Am.Vs.Sp))
grid.newpage()
grid.table(df_Am.Vs.Sp)
plotMD(et.Am.Vs.Sp)

# Volcano plot
volcanoData <- cbind(et.Am.Vs.Sp$table$logFC, -log10(et.Am.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Am.Vs.Sp$table[,"PValue"], breaks=50, main="Amygdala vs Spinal Cord p-value frequency histogram")

# Amygdala vs Substantia Nigra
et.Am.Vs.Sn <- exactTest(y, c("Substantia_Nigra", "Amygdala"))
df_Am.Vs.Sn <- summary(decideTests(et.Am.Vs.Sn))
grid.newpage()
grid.table(df_Am.Vs.Sn)
plotMD(et.Am.Vs.Sn)

# Volcano plot
volcanoData <- cbind(et.Am.Vs.Sn$table$logFC, -log10(et.Am.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Am.Vs.Sn$table[,"PValue"], breaks=50, main="Amygdala vs Substantia_Nigra p-value frequency histogram")

# Anterior vs Caudate
et.At.Vs.Ca <- exactTest(y, c("Caudate", "Anterior"))
df_At.Vs.Ca <- summary(decideTests(et.At.Vs.Ca))
grid.newpage()
grid.table(df_At.Vs.Ca)
plotMD(et.At.Vs.Ca)

# Volcano plot
volcanoData <- cbind(et.At.Vs.Ca$table$logFC, -log10(et.At.Vs.Ca$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Caudate')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.At.Vs.Ca$table[,"PValue"], breaks=50, main="Anterior vs Caudate p-value frequency histogram")

# Anterior vs Cerebellar
et.At.Vs.Ce <- exactTest(y, c("Cerebellar", "Anterior"))
df_At.Vs.Ce <- summary(decideTests(et.At.Vs.Ce))
grid.newpage()
grid.table(df_At.Vs.Ce)
plotMD(et.At.Vs.Ce)

# Volcano plot
volcanoData <- cbind(et.At.Vs.Ce$table$logFC, -log10(et.At.Vs.Ce$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Cerebellar')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.At.Vs.Ce$table[,"PValue"], breaks=50, main="Anterior vs Cerebellar p-value frequency histogram")

# Anterior vs Cerebellum
et.At.Vs.Cm <- exactTest(y, c("Cerebellum", "Anterior"))
df_At.Vs.Cm <- summary(decideTests(et.At.Vs.Cm))
grid.newpage()
grid.table(df_At.Vs.Cm)
plotMD(et.At.Vs.Cm)

# Volcano plot
volcanoData <- cbind(et.At.Vs.Cm$table$logFC, -log10(et.At.Vs.Cm$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Cerebellum')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.At.Vs.Cm$table[,"PValue"], breaks=50, main="Anterior vs Cerebellum p-value frequency histogram")

# Anterior vs Cortex
et.At.Vs.Co <- exactTest(y, c("Cortex", "Anterior"))
df_At.Vs.Co <- summary(decideTests(et.At.Vs.Co))
grid.newpage()
grid.table(df_At.Vs.Co)
plotMD(et.At.Vs.Co)

# Volcano plot
volcanoData <- cbind(et.At.Vs.Co$table$logFC, -log10(et.At.Vs.Co$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.At.Vs.Co$table[,"PValue"], breaks=50, main="Anterior vs Cortex p-value frequency histogram")

# Anterior vs Frontal Cortex
et.At.Vs.Fc <- exactTest(y, c("Frontal_Cortex", "Anterior"))
df_At.Vs.Fc <- summary(decideTests(et.At.Vs.Fc))
grid.newpage()
grid.table(df_At.Vs.Fc)
plotMD(et.At.Vs.Fc)

# Volcano plot
volcanoData <- cbind(et.At.Vs.Fc$table$logFC, -log10(et.At.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Frontal_Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.At.Vs.Fc$table[,"PValue"], breaks=50, main="Anterior vs Frontal Cortex p-value frequency histogram")

# Anterior vs Hippocampus
et.At.Vs.Hp <- exactTest(y, c("Hippocampus", "Anterior"))
df_At.Vs.Hp <- summary(decideTests(et.At.Vs.Hp))
grid.newpage()
grid.table(df_At.Vs.Hp)
plotMD(et.At.Vs.Hp)

# Volcano plot
volcanoData <- cbind(et.At.Vs.Hp$table$logFC, -log10(et.At.Vs.Hp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.At.Vs.Hp$table[,"PValue"], breaks=50, main="Anterior vs Hippocampus p-value frequency histogram")

# Anterior vs Hypothalamus
et.At.Vs.Hy <- exactTest(y, c("Hypothalamus", "Anterior"))
df_At.Vs.Hy <- summary(decideTests(et.At.Vs.Hy))
grid.newpage()
grid.table(df_At.Vs.Hy)
plotMD(et.At.Vs.Hy)

# Volcano plot
volcanoData <- cbind(et.At.Vs.Hy$table$logFC, -log10(et.At.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.At.Vs.Hp$table[,"PValue"], breaks=50, main="Anterior vs Hypothalamus p-value frequency histogram")

# Anterior vs Nucleus Accumbens
et.At.Vs.Nc <- exactTest(y, c("Nucleus_Accumbens", "Anterior"))
df_At.Vs.Nc <- summary(decideTests(et.At.Vs.Nc))
grid.newpage()
grid.table(df_At.Vs.Nc)
plotMD(et.At.Vs.Nc)

# Volcano plot
volcanoData <- cbind(et.At.Vs.Nc$table$logFC, -log10(et.At.Vs.Nc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.At.Vs.Nc$table[,"PValue"], breaks=50, main="Anterior vs Nucleus Accumbens p-value frequency histogram")

# Anterior vs Putamen
et.At.Vs.Pu <- exactTest(y, c("Putamen", "Anterior"))
df_At.Vs.Pu <- summary(decideTests(et.At.Vs.Pu))
grid.newpage()
grid.table(df_At.Vs.Pu)
plotMD(et.At.Vs.Pu)

# Volcano plot
volcanoData <- cbind(et.At.Vs.Pu$table$logFC, -log10(et.At.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.At.Vs.Pu$table[,"PValue"], breaks=50, main="Anterior vs Putamen p-value frequency histogram")

# Anterior vs Spinal Cord
et.At.Vs.Sp <- exactTest(y, c("Spinal_Cord", "Anterior"))
df_At.Vs.Sp <- summary(decideTests(et.At.Vs.Sp))
grid.newpage()
grid.table(df_At.Vs.Sp)
plotMD(et.At.Vs.Sp)

# Volcano plot
volcanoData <- cbind(et.At.Vs.Sp$table$logFC, -log10(et.At.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.At.Vs.Sp$table[,"PValue"], breaks=50, main="Anterior vs Spinal Cord p-value frequency histogram")

# Anterior vs Substantia Nigra
et.At.Vs.Sn <- exactTest(y, c("Substantia_Nigra", "Anterior"))
df_At.Vs.Sn <- summary(decideTests(et.At.Vs.Sn))
grid.newpage()
grid.table(df_At.Vs.Sn)
plotMD(et.At.Vs.Sn)

# Volcano plot
volcanoData <- cbind(et.At.Vs.Sn$table$logFC, -log10(et.At.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.At.Vs.Sn$table[,"PValue"], breaks=50, main="Anterior vs Substantia Nigra p-value frequency histogram")

# Caudate vs Cerebellar
et.Ca.Vs.Ce <- exactTest(y, c("Cerebellar", "Caudate"))
df_Ca.Vs.Ce <- summary(decideTests(et.Ca.Vs.Ce))
grid.newpage()
grid.table(df_Ca.Vs.Ce)
plotMD(et.Ca.Vs.Ce)

# Volcano plot
volcanoData <- cbind(et.Ca.Vs.Ce$table$logFC, -log10(et.Ca.Vs.Ce$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Cerebellar')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ca.Vs.Ce$table[,"PValue"], breaks=50, main="Caudate vs Cerebellar p-value frequency histogram")

# Caudate vs Cerebellum
et.Ca.Vs.Cm <- exactTest(y, c("Cerebellum", "Caudate"))
df_Ca.Vs.Cm <- summary(decideTests(et.Ca.Vs.Cm))
grid.newpage()
grid.table(df_Ca.Vs.Cm)
plotMD(et.At.Vs.Sn)

# Volcano plot
volcanoData <- cbind(et.Ca.Vs.Cm$table$logFC, -log10(et.Ca.Vs.Cm$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Cerebellum')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ca.Vs.Cm$table[,"PValue"], breaks=50, main="Caudate vs Cerebellum p-value frequency histogram")

# Caudate vs Cortex
et.Ca.Vs.Co <- exactTest(y, c("Cortex", "Caudate"))
df_Ca.Vs.Co <- summary(decideTests(et.Ca.Vs.Co))
grid.newpage()
grid.table(df_Ca.Vs.Co)
plotMD(et.Ca.Vs.Co)

# Volcano plot
volcanoData <- cbind(et.Ca.Vs.Co$table$logFC, -log10(et.Ca.Vs.Co$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ca.Vs.Co$table[,"PValue"], breaks=50, main="Caudate vs Cortex p-value frequency histogram")

# Caudate vs Frontal Cortex
et.Ca.Vs.Fc <- exactTest(y, c("Frontal_Cortex", "Caudate"))
df_Ca.Vs.Fc <- summary(decideTests(et.Ca.Vs.Fc))
grid.newpage()
grid.table(df_Ca.Vs.Fc)
plotMD(et.Ca.Vs.Fc)

# Volcano plot
volcanoData <- cbind(et.Ca.Vs.Fc$table$logFC, -log10(et.Ca.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Frontal_Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ca.Vs.Fc$table[,"PValue"], breaks=50, main="Caudate vs Frontal Cortex p-value frequency histogram")

# Caudate vs Hipocampus
et.Ca.Vs.Hp <- exactTest(y, c("Hippocampus", "Caudate"))
df_Ca.Vs.Hp <- summary(decideTests(et.Ca.Vs.Hp))
grid.newpage()
grid.table(df_Ca.Vs.Hp)
plotMD(et.Ca.Vs.Hp)

# Volcano plot
volcanoData <- cbind(et.Ca.Vs.Hp$table$logFC, -log10(et.Ca.Vs.Hp[$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ca.Vs.Hp$table[,"PValue"], breaks=50, main="Caudate vs Hippocampus p-value frequency histogram")

# Caudate vs Hypothalamus
et.Ca.Vs.Hy <- exactTest(y, c("Hypothalamus", "Caudate"))
df_Ca.Vs.Hy <- summary(decideTests(et.Ca.Vs.Hy))
grid.newpage()
grid.table(df_Ca.Vs.Hy)
plotMD(et.Ca.Vs.Hy)

# Volcano plot
volcanoData <- cbind(et.Ca.Vs.Hy$table$logFC, -log10(et.Ca.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ca.Vs.Hy$table[,"PValue"], breaks=50, main="Caudate vs Hypothalamus p-value frequency histogram")

# Caudate vs Nucleus Accumbens
et.Ca.Vs.Na <- exactTest(y, c("Nucleus_Accumebns", "Caudate"))
df_Ca.Vs.Na <- summary(decideTests(et.Ca.Vs.Na))
grid.newpage()
grid.table(df_Ca.Vs.Na)
plotMD(et.Ca.Vs.Na)

# Volcano plot
volcanoData <- cbind(et.Ca.Vs.Na$table$logFC, -log10(et.Ca.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ca.Vs.Na$table[,"PValue"], breaks=50, main="Caudate vs Nucleus Accumbens p-value frequency histogram")

# Caudate vs Putamen
et.Ca.Vs.Pu <- exactTest(y, c("Putamen", "Caudate"))
df_Ca.Vs.Pu <- summary(decideTests(et.Ca.Vs.Pu))
grid.newpage()
grid.table(df_Ca.Vs.Pu)
plotMD(et.Ca.Vs.Pu)

# Volcano plot
volcanoData <- cbind(et.Ca.Vs.Pu$table$logFC, -log10(et.Ca.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ca.Vs.Pu$table[,"PValue"], breaks=50, main="Caudate vs Putamen p-value frequency histogram")

# Caudate vs Spinal Cord
et.Ca.Vs.Sp <- exactTest(y, c("Spinal_Cord", "Caudate"))
df_Ca.Vs.Sp <- summary(decideTests(et.Ca.Vs.Sp))
grid.newpage()
grid.table(df_Ca.Vs.Sp)
plotMD(et.Ca.Vs.Sp)

# Volcano plot
volcanoData <- cbind(et.Ca.Vs.Sp$table$logFC, -log10(et.Ca.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ca.Vs.Sp$table[,"PValue"], breaks=50, main="Caudate vs Spinal Cord p-value frequency histogram")

# Caudate vs Substantia Nigra
et.Ca.Vs.Sn <- exactTest(y, c("Substantia_Nigra", "Caudate"))
df_Ca.Vs.Sn <- summary(decideTests(et.Ca.Vs.Sn))
grid.newpage()
grid.table(df_Ca.Vs.Sn)
plotMD(et.Ca.Vs.Sn)

# Volcano plot
volcanoData <- cbind(et.Ca.Vs.Sn$table$logFC, -log10(et.Ca.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ca.Vs.Sn$table[,"PValue"], breaks=50, main="Caudate vs Substantia Nigra p-value frequency histogram")

# Cerebellar vs Cerebellum
et.Ce.Vs.Cm <- exactTest(y, c("Cerebellum", "Cerebellar"))
df_Ce.Vs.Cm <- summary(decideTests(et.Ce.Vs.Cm))
grid.newpage()
grid.table(df_Ce.Vs.Cm)
plotMD(et.Ce.Vs.Cm)

# Volcano plot
volcanoData <- cbind(et.Ce.Vs.Cm$table$logFC, -log10(et.Ce.Vs.Cm$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Cerebellum')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ce.Vs.Cm$table[,"PValue"], breaks=50, main="Cerebellar vs Cerebellum p-value frequency histogram")

# Cerebellar vs Cortex
et.Ce.Vs.Co <- exactTest(y, c("Cortex", "Cerebellar"))
df_Ce.Vs.Co <- summary(decideTests(et.Ce.Vs.Co))
grid.newpage()
grid.table(df_Ce.Vs.Co)
plotMD(et.Ce.Vs.Co)

# Volcano plot
volcanoData <- cbind(et.Ce.Vs.Co$table$logFC, -log10(et.Ce.Vs.Co$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ce.Vs.Co$table[,"PValue"], breaks=50, main="Cerebellar vs Cortex p-value frequency histogram")

# Cerebellar vs Frontal Cortex
et.Ce.Vs.Fc <- exactTest(y, c("Frontal_Cortex", "Cerebellar"))
df_Ce.Vs.Fc <- summary(decideTests(et.Ce.Vs.Fc))
grid.newpage()
grid.table(df_Ce.Vs.Fc)
plotMD(et.Ce.Vs.Fc)

# Volcano plot
volcanoData <- cbind(et.Ce.Vs.Fc$table$logFC, -log10(et.Ce.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Frontal_Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ce.Vs.Fc$table[,"PValue"], breaks=50, main="Cerebellar vs Frontal Cortex p-value frequency histogram")

# Cerebellar vs Hippocampus
et.Ce.Vs.Hp <- exactTest(y, c("Hippocampus", "Cerebellar"))
df_Ce.Vs.Hp <- summary(decideTests(et.Ce.Vs.Hp))
grid.newpage()
grid.table(df_Ce.Vs.Hp)
plotMD(et.Ce.Vs.Hp)

# Volcano plot
volcanoData <- cbind(et.Ce.Vs.Hp$table$logFC, -log10(et.Ce.Vs.Hp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ce.Vs.Hp$table[,"PValue"], breaks=50, main="Cerebellar vs Hippocampus p-value frequency histogram")

# Cerebellar vs Hypothalamus
et.Ce.Vs.Hy <- exactTest(y, c("Hypothalamus", "Cerebellar"))
df_Ce.Vs.Hy <- summary(decideTests(et.Ce.Vs.Hy))
grid.newpage()
grid.table(df_Ce.Vs.Hy)
plotMD(et.Ce.Vs.Hy)

# Volcano plot
volcanoData <- cbind(et.Ce.Vs.Hy$table$logFC, -log10(et.Ce.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ce.Vs.Hy$table[,"PValue"], breaks=50, main="Cerebellar vs Hypothalamus p-value frequency histogram")

# Cerebellar vs Nucleus Accumbens
et.Ce.Vs.Na <- exactTest(y, c("Nucleus_Accumbens", "Cerebellar"))
df_Ce.Vs.Na <- summary(decideTests(et.Ce.Vs.Na))
grid.newpage()
grid.table(df_Ce.Vs.Na)
plotMD(et.Ce.Vs.Na)

# Volcano plot
volcanoData <- cbind(et.Ce.Vs.Na$table$logFC, -log10(et.Ce.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ce.Vs.Na$table[,"PValue"], breaks=50, main="Cerebellar vs Nucleus Accumbens p-value frequency histogram")

# Cerebellar vs Putamen
et.Ce.Vs.Pu <- exactTest(y, c("Putamen", "Cerebellar"))
df_Ce.Vs.Pu <- summary(decideTests(et.Ce.Vs.Pu))
grid.newpage()
grid.table(df_Ce.Vs.Pu)
plotMD(et.Ce.Vs.Pu)

# Volcano plot
volcanoData <- cbind(et.Ce.Vs.Pu$table$logFC, -log10(et.Ce.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ce.Vs.Pu$table[,"PValue"], breaks=50, main="Cerebellar vs Putamen p-value frequency histogram")

# Cerebellar vs Spinal Cord
et.Ce.Vs.Sp <- exactTest(y, c("Spinal_Cord", "Cerebellar"))
df_Ce.Vs.Sp <- summary(decideTests(et.Ce.Vs.Sp))
grid.newpage()
grid.table(df_Ce.Vs.Sp)
plotMD(et.Ce.Vs.Sp)

# Volcano plot
volcanoData <- cbind(et.Ce.Vs.Sp$table$logFC, -log10(et.Ce.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ce.Vs.Sp$table[,"PValue"], breaks=50, main="Cerebellar vs Spinal Cord p-value frequency histogram")

# Cerebellar vs Substantia Nigra
et.Ce.Vs.Sn <- exactTest(y, c("Substantia_Nigra", "Cerebellar"))
df_Ce.Vs.Sn <- summary(decideTests(et.Ce.Vs.Sn))
grid.newpage()
grid.table(df_Ce.Vs.Sn)
plotMD(et.Ce.Vs.Sn)

# Volcano plot
volcanoData <- cbind(et.Ce.Vs.Sn$table$logFC, -log10(et.Ce.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Ce.Vs.Sn$table[,"PValue"], breaks=50, main="Cerebellar vs Substantia Nigra p-value frequency histogram")

# Cortex vs Frontal Cortex
et.Co.Vs.Fc <- exactTest(y, c("Frontal_Cortex", "Cortex"))
df_Co.Vs.Fc <- summary(decideTests(et.Co.Vs.Fc))
grid.newpage()
grid.table(df_Co.Vs.Fc)
plotMD(et.Co.Vs.Fc)

# Volcano plot
volcanoData <- cbind(et.Co.Vs.Fc$table$logFC, -log10(et.Co.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortexr-1*Frontal_Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Co.Vs.Fc$table[,"PValue"], breaks=50, main="Cortex vs Frontal Cortex p-value frequency histogram")

# Cortex vs Hippocampus
et.Co.Vs.Hp <- exactTest(y, c("Hippocampus", "Cortex"))
df_Co.Vs.Hp <- summary(decideTests(et.Co.Vs.Hp))
grid.newpage()
grid.table(df_Co.Vs.Hp)
plotMD(et.Co.Vs.Hp)

# Volcano plot
volcanoData <- cbind(et.Co.Vs.Hp$table$logFC, -log10(et.Co.Vs.Hp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Co.Vs.Hp$table[,"PValue"], breaks=50, main="Cortex vs Hippocampus p-value frequency histogram")

# Cortex vs Hypothalamus
et.Co.Vs.Hy <- exactTest(y, c("Hypothalamus", "Cortex"))
df_Co.Vs.Hy <- summary(decideTests(et.Co.Vs.Hy))
grid.newpage()
grid.table(df_Co.Vs.Hy)
plotMD(et.Co.Vs.Hy)

# Volcano plot
volcanoData <- cbind(et.Co.Vs.Fc$table$logFC, -log10(et.Co.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Co.Vs.Hy$table[,"PValue"], breaks=50, main="Cortex vs Hypothalamus p-value frequency histogram")

# Cortex vs Nucleus Accumbens
et.Co.Vs.Na <- exactTest(y, c("Nucleus_Accumbens", "Cortex"))
df_Co.Vs.Na <- summary(decideTests(et.Co.Vs.Na))
grid.newpage()
grid.table(df_Co.Vs.Na)
plotMD(et.Co.Vs.Na)

# Volcano plot
volcanoData <- cbind(et.Co.Vs.Na$table$logFC, -log10(et.Co.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Co.Vs.Na$table[,"PValue"], breaks=50, main="Cortex vs Nucleus Accumbens p-value frequency histogram")

# Cortex vs Putamen
et.Co.Vs.Pu <- exactTest(y, c("Putamen", "Cortex"))
df_Co.Vs.Pu <- summary(decideTests(et.Co.Vs.Pu))
grid.newpage()
grid.table(df_Co.Vs.Pu)
plotMD(et.Co.Vs.Pu)

# Volcano plot
volcanoData <- cbind(et.Co.Vs.Pu$table$logFC, -log10(et.Co.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Co.Vs.Pu$table[,"PValue"], breaks=50, main="Cortex vs Putamen p-value frequency histogram")

# Cortex vs Spinal Cord
et.Co.Vs.Sp <- exactTest(y, c("Spinal_Cord", "Cortex"))
df_Co.Vs.Sp <- summary(decideTests(et.Co.Vs.Sp))
grid.newpage()
grid.table(df_Co.Vs.Sp)
plotMD(et.Co.Vs.Sp)

# Volcano plot
volcanoData <- cbind(et.Co.Vs.Sp$table$logFC, -log10(et.Co.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Co.Vs.Sp$table[,"PValue"], breaks=50, main="Cortex vs Spinal Cord p-value frequency histogram")

# Cortex vs Substantia Nigra
et.Co.Vs.Sn <- exactTest(y, c("Substantia_Nigra", "Cortex"))
df_Co.Vs.Sn <- summary(decideTests(et.Co.Vs.Sn))
grid.newpage()
grid.table(df_Co.Vs.Sn)
plotMD(et.Co.Vs.Sn)

# Volcano plot
volcanoData <- cbind(et.Co.Vs.Sn$table$logFC, -log10(et.Co.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Co.Vs.Sn$table[,"PValue"], breaks=50, main="Cortex vs Substantia Nigra p-value frequency histogram")

# Frontal Cortex vs Hippocampus
et.Fc.Vs.Hp <- exactTest(y, c("Hippocampus", "Frontal_Cortex"))
df_Fc.Vs.Hp <- summary(decideTests(et.Fc.Vs.Hp))
grid.newpage()
grid.table(df_Fc.Vs.Hp)
plotMD(et.Fc.Vs.Hp)

# Volcano plot
volcanoData <- cbind(et.Fc.Vs.Hp$table$logFC, -log10(et.Fc.Vs.Hp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Fc.Vs.Hp$table[,"PValue"], breaks=50, main="Frontal Cortex vs Hippocampus p-value frequency histogram")

# Frontal Cortex vs Hypothalamus 
et.Fc.Vs.Hy <- exactTest(y, c("Hypothalamus", "Frontal_Cortex"))
df_Fc.Vs.Hy <- summary(decideTests(et.Fc.Vs.Hy))
grid.newpage()
grid.table(df_Fc.Vs.Hy)
plotMD(et.Fc.Vs.Hy)

# Volcano plot
volcanoData <- cbind(et.Fc.Vs.Hy$table$logFC, -log10(et.Fc.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Fc.Vs.Hy$table[,"PValue"], breaks=50, main="Frontal Cortex vs Hypothalamus p-value frequency histogram")

# Frontal Cortex vs Nucleus Accumbens
et.Fc.Vs.Na <- exactTest(y, c("Nucleus_Accumbens", "Frontal_Cortex"))
df_Fc.Vs.Na <- summary(decideTests(et.Fc.Vs.Na))
grid.newpage()
grid.table(df_Fc.Vs.Na)
plotMD(et.Fc.Vs.Na)

# Volcano plot
volcanoData <- cbind(et.Fc.Vs.Na$table$logFC, -log10(et.Fc.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Fc.Vs.Na$table[,"PValue"], breaks=50, main="Frontal Cortex vs Nucleus Accumbens p-value frequency histogram")

# Frontal Cortex vs Putamen
et.Fc.Vs.Pu <- exactTest(y, c("Putamen", "Frontal_Cortex"))
df_Fc.Vs.Pu <- summary(decideTests(et.Fc.Vs.Pu))
grid.newpage()
grid.table(df_Fc.Vs.Pu)
plotMD(et.Fc.Vs.Pu)

# Volcano plot
volcanoData <- cbind(et.Fc.Vs.Pu$table$logFC, -log10(et.Fc.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Fc.Vs.Pu$table[,"PValue"], breaks=50, main="Frontal Cortex vs Putamen p-value frequency histogram")

# Frontal Corte vs Spinal Cord
et.Fc.Vs.Sp <- exactTest(y, c("Spinal_Cord", "Frontal_Cortex"))
df_Fc.Vs.Sp <- summary(decideTests(et.Fc.Vs.Sp))
grid.newpage()
grid.table(df_Fc.Vs.Sp)
plotMD(et.Fc.Vs.Sp)

# Volcano plot
volcanoData <- cbind(et.Fc.Vs.Sp$table$logFC, -log10(et.Fc.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Fc.Vs.Sp$table[,"PValue"], breaks=50, main="Frontal Cortex vs Spinal Cord p-value frequency histogram")

# Frontal Cortex vs Substantia Nigra
et.Fc.Vs.Sn <- exactTest(y, c("Substantia_Nigra", "Frontal_Cortex"))
df_Fc.Vs.Sn <- summary(decideTests(et.Fc.Vs.Sn))
grid.newpage()
grid.table(df_Fc.Vs.Sn)
plotMD(et.Fc.Vs.Sn)

# Volcano plot
volcanoData <- cbind(et.Fc.Vs.Sn$table$logFC, -log10(et.Fc.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Fc.Vs.Sn$table[,"PValue"], breaks=50, main="Frontal Cortex vs Substantia Nigra p-value frequency histogram")

# Hippocampus vs Hypothalamus
et.Hp.Vs.Hy <- exactTest(y, c("Hypothalamus", "Hippocampus"))
df_Hp.Vs.Hy <- summary(decideTests(et.Hp.Vs.Hy))
grid.newpage()
grid.table(df_Hp.Vs.Hy)
plotMD(et.Hp.Vs.Hy)

# Volcano plot
volcanoData <- cbind(et.Hp.Vs.Hy$table$logFC, -log10(et.Hp.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Hp.Vs.Hy$table[,"PValue"], breaks=50, main="Hippocampus vs Hypothalamus p-value frequency histogram")

# Hippocampus vs Nucleus Accumbens
et.Hp.Vs.Na <- exactTest(y, c("Nucleus_Accumbens", "Hippocampus"))
df_Hp.Vs.Na <- summary(decideTests(et.Hp.Vs.Na))
grid.newpage()
grid.table(df_Hp.Vs.Na)
plotMD(et.Hp.Vs.Na)

# Volcano plot
volcanoData <- cbind(et.Hp.Vs.Na$table$logFC, -log10(et.Hp.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Hp.Vs.Na$table[,"PValue"], breaks=50, main="Hippocampus vs Nucleus Accumbens p-value frequency histogram")

# Hippocampus vs Putamen
et.Hp.Vs.Pu <- exactTest(y, c("Putamen", "Hippocampus"))
df_Hp.Vs.Pu <- summary(decideTests(et.Hp.Vs.Pu))
grid.newpage()
grid.table(df_Hp.Vs.Pu)
plotMD(et.Hp.Vs.Pu)

# Volcano plot
volcanoData <- cbind(et.Hp.Vs.Pu$table$logFC, -log10(et.Hp.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Hp.Vs.Pu$table[,"PValue"], breaks=50, main="Hippocampus vs Putamen p-value frequency histogram")

# Hippocampus vs Spinal Cord
et.Hp.Vs.Sp <- exactTest(y, c("Spinal_Cord", "Hippocampus"))
df_Hp.Vs.Sp <- summary(decideTests(et.Hp.Vs.Sp))
grid.newpage()
grid.table(df_Hp.Vs.Sp)
plotMD(et.Hp.Vs.Sp)

# Volcano plot
volcanoData <- cbind(et.Hp.Vs.Sp$table$logFC, -log10(et.Hp.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Hp.Vs.Sp$table[,"PValue"], breaks=50, main="Hippocampus vs Spinal Cord p-value frequency histogram")

# Hippocampus vs Substantia Nigra 
et.Hp.Vs.Sn <- exactTest(y, c("Substantia_Nigra", "Hippocampus"))
df_Hp.Vs.Sn <- summary(decideTests(et.Hp.Vs.Sn))
grid.newpage()
grid.table(df_Hp.Vs.Sn)
plotMD(et.Hp.Vs.Sn)

# Volcano plot
volcanoData <- cbind(et.Hp.Vs.Sn$table$logFC, -log10(et.Hp.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Hp.Vs.Sn$table[,"PValue"], breaks=50, main="Hippocampus vs Substantia Nigra p-value frequency histogram")

# Hypothalamus vs Nucleus Accumbens
et.Hy.Vs.Na <- exactTest(y, c("Nucleus_Accumbens", "Hypothalamus"))
df_Hy.Vs.Na <- summary(decideTests(et.Hy.Vs.Na))
grid.newpage()
grid.table(df_Hy.Vs.Na)
plotMD(et.Hy.Vs.Na)

# Volcano plot
volcanoData <- cbind(et.Hy.Vs.Na$table$logFC, -log10(et.Hy.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hypothalamus-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Hy.Vs.Na$table[,"PValue"], breaks=50, main="Hypothalamus vs Nucleus Accumbens p-value frequency histogram")

# Hypothalamus vs Putamen
et.Hy.Vs.Pu <- exactTest(y, c("Putamen", "Hypothalamus"))
df_Hy.Vs.Pu <- summary(decideTests(et.Hy.Vs.Pu))
grid.newpage()
grid.table(df_Hy.Vs.Pu)
plotMD(et.Hy.Vs.Pu)

# Volcano plot
volcanoData <- cbind(et.Hy.Vs.Pu$table$logFC, -log10(et.Hy.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hypothalamus-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Hy.Vs.Pu$table[,"PValue"], breaks=50, main="Hypothalamus vs Putamen p-value frequency histogram")

# Hypothalamus vs Spinal Cord
et.Hy.Vs.Sp <- exactTest(y, c("Spinal_Cord", "Hypothalamus"))
df_Hy.Vs.Sp <- summary(decideTests(et.Hy.Vs.Sp))
grid.newpage()
grid.table(df_Hy.Vs.Sp)
plotMD(et.Hy.Vs.Sp)

# Volcano plot
volcanoData <- cbind(et.Hy.Vs.Sp$table$logFC, -log10(et.Hy.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hypothalamus-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Hy.Vs.Sp$table[,"PValue"], breaks=50, main="Hypothalamus vs Spinal Cord p-value frequency histogram")

# Hypothalamus Substantia Nigra
et.Hy.Vs.Sn <- exactTest(y, c("Substantia_Nigra", "Hypothalamus"))
df_Hy.Vs.Sn <- summary(decideTests(et.Hy.Vs.Sn))
grid.newpage()
grid.table(df_Hy.Vs.Sn)
plotMD(et.Hy.Vs.Sn)

# Volcano plot
volcanoData <- cbind(et.Hy.Vs.Sn$table$logFC, -log10(et.Hy.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hypothalamus-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Hy.Vs.Sn$table[,"PValue"], breaks=50, main="Hypothalamus vs Substantia Nigra p-value frequency histogram")

# Nucleus Accumbens vs Putamen
et.Na.Vs.Pu <- exactTest(y, c("Putamen", "Nucleus_Accumbens"))
df_Na.Vs.Pu <- summary(decideTests(et.Na.Vs.Pu))
grid.newpage()
grid.table(df_Na.Vs.Pu)
plotMD(et.Na.Vs.Pu)

# Volcano plot
volcanoData <- cbind(et.Na.Vs.Pu$table$logFC, -log10(et.Na.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Nucleus_Accumbens-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Na.Vs.Pu$table[,"PValue"], breaks=50, main="Nucleus Accumbens vs Putamen p-value frequency histogram")

# Nucleus Accumbens vs Spinal Cord
et.Na.Vs.Sp <- exactTest(y, c("Spinal_Cord", "Nucleus_Accumbens"))
df_Na.Vs.Sp <- summary(decideTests(et.Na.Vs.Sp))
grid.newpage()
grid.table(df_Na.Vs.Sp)
plotMD(et.Na.Vs.Sp)

# Volcano plot
volcanoData <- cbind(et.Na.Vs.Sp$table$logFC, -log10(et.Na.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Nucleus_Accumbens-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Na.Vs.Sp$table[,"PValue"], breaks=50, main="Nucleus Accumbens vs Spinal Cord p-value frequency histogram")

# Nucleus Accumbens vs Substantia Nigra
et.Na.Vs.Sn <- exactTest(y, c("Substantia_Nigra", "Nucleus_Accumbens"))
df_Na.Vs.Sn <- summary(decideTests(et.Na.Vs.Sn))
grid.newpage()
grid.table(df_na.Vs.Sn)
plotMD(et.Na.Vs.Sn)

# Volcano plot
volcanoData <- cbind(et.Na.Vs.Sn$table$logFC, -log10(et.Na.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Nucleus_Accumbens-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Na.Vs.Sn$table[,"PValue"], breaks=50, main="Nucleus Accumbens vs Substantia Nigra p-value frequency histogram")

# Putamen vs Spinal Cord
et.Pu.Vs.Sp <- exactTest(y, c("Spinal_Cord", "Putamen"))
df_Pu.Vs.Sp <- summary(decideTests(et.Pu.Vs.Sp))
grid.newpage()
grid.table(df_Pu.Vs.Sp)
plotMD(et.Pu.Vs.Sp)

# Volcano plot
volcanoData <- cbind(et.Pu.Vs.Sp$table$logFC, -log10(et.Pu.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Putamen-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Pu.Vs.Sp$table[,"PValue"], breaks=50, main="Putamen vs Spinal Cord p-value frequency histogram")

# Putamen vs Substantia Nigra
et.Pu.Vs.Sn <- exactTest(y, c("Substantia_Nigra", "Putamen"))
df_Pu.Vs.Sn <- summary(decideTests(et.Pu.Vs.Sn))
grid.newpage()
grid.table(df_Pu.Vs.Sn)
plotMD(et.Pu.Vs.Sn)

# Volcano plot
volcanoData <- cbind(et.Pu.Vs.Sn$table$logFC, -log10(et.Pu.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Putamen-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Pu.Vs.Sn$table[,"PValue"], breaks=50, main="Putamen vs Substantia Nigra p-value frequency histogram")

# Spinal Cord vs Substantia Nigra
et.Sp.Vs.Sn <- exactTest(y, c("Substantia_Nigra", "Spinal_Cord"))
df_Sp.Vs.Sn <- summary(decideTests(et.Sp.Vs.Sn))
grid.newpage()
grid.table(df_Sp.Vs.Sn)
plotMD(et.Sp.Vs.Sn)

# Volcano plot
volcanoData <- cbind(et.Sp.Vs.Sn$table$logFC, -log10(et.Sp.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Spinal_Cord-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(et.Sp.Vs.Sn$table[,"PValue"], breaks=50, main="Spinal Cord vs Substantia Nigra p-value frequency histogram")

dev.off()

