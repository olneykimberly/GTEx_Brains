# This script looks at differential gene expression between each brain tissue type.

METADATA <- '/scratch/mjpete11/GTEx/Metadata.csv'
COUNT_MATRIX <- '/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrix.tsv'
PLOT_DIR =  '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Tissue/Gene/Normalized/GLM_ratio_Test/Plots/'
FILE_NAME = 'Unnormalized_GLM_Ratio_Test.pdf'

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

# Test for DGX using liklihood ratio test.
fit <- glmFit(y, design, robust=TRUE)

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


# Test for DGX using quasi-liklihood method.

# Plot
setwd(PLOT_DIR)

pdf(FILE_NAME)

# Amygdala vs Anterior
lrt.Am.Vs.At <- glmLRT(fit, contrast=my.contrasts[,"Am.Vs.At"])
df_Am.Vs.At <- summary(decideTests(lrt.Am.Vs.At))
grid.newpage()
grid.table(df_Am.Vs.At)
plotMD(lrt.Am.Vs.At)

# Volcano plot
volcanoData <- cbind(lrt.Am.Vs.At$table$logFC, -log10(lrt.Am.Vs.At$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Anterior')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Am.Vs.At$table[,"PValue"], breaks=50, main="Amygdala vs Anterior p-value frequency histogram")

# Amygdala vs Caduate
lrt.Am.Vs.Ca <- glmLRT(fit, contrast=my.contrast[,"Am.Vs.Ca"])
df_Am.Vs.Ca <- summary(decideTests(lrt.Am.Vs.Ca))
grid.newpage()
grid.table(df_Am.Vs.Ca)
plotMD(lrt.Am.Vs.Ca)

# Volcano plot
volcanoData <- cbind(lrt.Am.Vs.Ca$table$logFC, -log10(lrt.Am.Vs.Ca$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Caudate')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Am.Vs.Ca$table[,"PValue"], breaks=50, main="Amygdala vs Caudate p-value frequency histogram")

# Amygdala vs Cerebellar
lrt.Am.Vs.Ce <- glmLRT(fit, contrast=my.contrasts[,"Am.Vs.Ce"])
df_Am.Vs.Ce <- summary(decideTests(lrt.Am.Vs.Ce))
grid.newpage()
grid.table(df_Am.Vs.Ce)
plotMD(lrt.Am.Vs.Ce)

# Volcano plot
volcanoData <- cbind(lrt.Am.Vs.Ce$table$logFC, -log10(lrt.Am.Vs.Ce$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Cerebellar')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Am.Vs.Ce$table[,"PValue"], breaks=50, main="Amygdala vs Cerebellar p-value frequency histogram")

# Amygdala vs Cerebellum
lrt.Am.Vs.Cm <- glmLRT(fit, contrast=my.contrasts[,"Am.Vs.Cm"])
df_Am.Vs.Cm <- summary(decideTests(lrt.Am.Vs.Cm))
grid.newpage()
grid.table(df_Am.Vs.Cm)
plotMD(lrt.Am.Vs.Cm)

# Volcano plot
volcanoData <- cbind(lrt.Am.Vs.Cm$table$logFC, -log10(lrt.Am.Vs.Cm$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Cerebellum')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Am.Vs.Cm$table[,"PValue"], breaks=50, main="Amygdala vs Cerebellum p-value frequency histogram")

# Amygdala vs Cortex
lrt.Am.Vs.Co <- glmLRT(fit, contrast=my.contrasts[,"Am.Vs.Co"])
df_Am.Vs.Co <- summary(decideTests(lrt.Am.Vs.Co))
grid.newpage()
grid.table(df_Am.Vs.Co)
plotMD(lrt.Am.Vs.Co)

# Volcano plot
volcanoData <- cbind(lrt.Am.Vs.Co$table$logFC, -log10(lrt.Am.Vs.Co$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Am.Vs.Co$table[,"PValue"], breaks=50, main="Amygdala vs Cortex p-value frequency histogram")

# Amygdala vs Frontal Cortex
lrt.Am.Vs.Fc <- glmLRT(fit, contrast=my.contrasts[,"Am.Vs.Fc"])
df_Am.Vs.Fc <- summary(decideTests(lrt.Am.Vs.Fc))
grid.newpage()
grid.table(df_Am.Vs.Fc)
plotMD(lrt.Am.Vs.Fc)

# Volcano plot
volcanoData <- cbind(lrt.Am.Vs.Fc$table$logFC, -log10(lrt.Am.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Frontal_Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Am.Vs.Fc$table[,"PValue"], breaks=50, main="Amygdala vs Frontal Cortex p-value frequency histogram")

# Amygdala vs Hippocampus
lrt.Am.Vs.Hp <- glmLRT(fit, contrast=my.contrasts[,"Am.Vs.Hp"])
summary(decideTests(lrt.Am.Vs.Hp))
plotMD(lrt.Am.Vs.Hp)

# Volcano plot
volcanoData <- cbind(lrt.Am.Vs.Hp$table$logFC, -log10(lrt.Am.Vs.Hp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Am.Vs.Hp$table[,"PValue"], breaks=50, main="Amygdala vs Hippocampus p-value frequency histogram")

# Amygdala vs Hypothalamus
lrt.Am.Vs.Hy <- glmLRT(fit, contrast=my.contrasts[,"Am.Vs.Hy"])
df_Am.Vs.Hy <- summary(decideTests(lrt.Am.Vs.Hy))
grid.newpage()
grid.table(df_Am.Vs.Hy)
plotMD(lrt.Am.Vs.Hy)

# Volcano plot
volcanoData <- cbind(lrt.Am.Vs.Hy$table$logFC, -log10(lrt.Am.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Am.Vs.Hy$table[,"PValue"], breaks=50, main="Amygdala vs Hypothalamus p-value frequency histogram")

# Amygdala vs Nucleus Accumbens
lrt.Am.Vs.Nc <- glmLRT(fit, contrast=my.contrasts[,"Am.Vs.Nc"])
df_Am.Vs.Nc <- summary(decideTests(lrt.Am.Vs.Nc))
grid.newpage()
grid.table(df_Am.Vs.Nc)
plotMD(lrt.Am.Vs.Nc)

# Volcano plot
volcanoData <- cbind(lrt.Am.Vs.Nc$table$logFC, -log10(lrt.Am.Vs.Nc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Am.Vs.Nc$table[,"PValue"], breaks=50, main="Amygdala vs Nucleus Accumbens p-value frequency histogram")

# Amygdala vs Putamen
lrt.Am.Vs.Pu <- glmLRT(fit, contrast=my.contrasts[,"Am.Vs.Pu"])
df_Am.Vs.Pu <- summary(decideTests(lrt.Am.Vs.Pu))
grid.newpage()
grid.table(df_Am.Vs.Pu)
plotMD(lrt.Am.Vs.Pu)

# Volcano plot
volcanoData <- cbind(lrt.Am.Vs.Pu$table$logFC, -log10(lrt.Am.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Am.Vs.Pu$table[,"PValue"], breaks=50, main="Amygdala vs Putamen p-value frequency histogram")

# Amygdala vs Spinal Cord
lrt.Am.Vs.Sp <- glmLRT(fit, contrast=my.contrasts[,"Am.Vs.Sp"])
df_Am.Vs.Sp <- summary(decideTests(lrt.Am.Vs.Sp))
grid.newpage()
grid.table(df_Am.Vs.Sp)
plotMD(lrt.Am.Vs.Sp)

# Volcano plot
volcanoData <- cbind(lrt.Am.Vs.Sp$table$logFC, -log10(lrt.Am.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Am.Vs.Sp$table[,"PValue"], breaks=50, main="Amygdala vs Spinal Cord p-value frequency histogram")

# Amygdala vs Substantia Nigra
lrt.Am.Vs.Sn <- glmLRT(fit, contrast=my.contrasts[,"Am.Vs.Sn"])
df_Am.Vs.Sn <- summary(decideTests(lrt.Am.Vs.Sn))
grid.newpage()
grid.table(df_Am.Vs.Sn)
plotMD(lrt.Am.Vs.Sn)

# Volcano plot
volcanoData <- cbind(lrt.Am.Vs.Sn$table$logFC, -log10(lrt.Am.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Am.Vs.Sn$table[,"PValue"], breaks=50, main="Amygdala vs Substantia_Nigra p-value frequency histogram")

# Anterior vs Caudate
lrt.At.Vs.Ca <- glmLRT(fit, contrast=my.contrasts[,"At.Vs.Ca"])
df_At.Vs.Ca <- summary(decideTests(lrt.At.Vs.Ca))
grid.newpage()
grid.table(df_At.Vs.Ca)
plotMD(lrt.At.Vs.Ca)

# Volcano plot
volcanoData <- cbind(lrt.At.Vs.Ca$table$logFC, -log10(lrt.At.Vs.Ca$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Caudate')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.At.Vs.Ca$table[,"PValue"], breaks=50, main="Anterior vs Caudate p-value frequency histogram")

# Anterior vs Cerebellar
lrt.At.Vs.Ce <- glmLRT(fit, contrast=my.contrasts[,"At.Vs.Ce"])
df_At.Vs.Ce <- summary(decideTests(lrt.At.Vs.Ce))
grid.newpage()
grid.table(df_At.Vs.Ce)
plotMD(lrt.At.Vs.Ce)

# Volcano plot
volcanoData <- cbind(lrt.At.Vs.Ce$table$logFC, -log10(lrt.At.Vs.Ce$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Cerebellar')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.At.Vs.Ce$table[,"PValue"], breaks=50, main="Anterior vs Cerebellar p-value frequency histogram")

# Anterior vs Cerebellum
lrt.At.Vs.Cm <- glmLRT(fit, contrast=my.contrasts[,"At.Vs.Cm"])
df_At.Vs.Cm <- summary(decideTests(lrt.At.Vs.Cm))
grid.newpage()
grid.table(df_At.Vs.Cm)
plotMD(lrt.At.Vs.Cm)

# Volcano plot
volcanoData <- cbind(lrt.At.Vs.Cm$table$logFC, -log10(lrt.At.Vs.Cm$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Cerebellum')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.At.Vs.Cm$table[,"PValue"], breaks=50, main="Anterior vs Cerebellum p-value frequency histogram")

# Anterior vs Cortex
lrt.At.Vs.Co <- glmLRT(fit, contrast=my.contrasts[,"At.Vs.Co"])
df_At.Vs.Co <- summary(decideTests(lrt.At.Vs.Co))
grid.newpage()
grid.table(df_At.Vs.Co)
plotMD(lrt.At.Vs.Co)

# Volcano plot
volcanoData <- cbind(lrt.At.Vs.Co$table$logFC, -log10(lrt.At.Vs.Co$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.At.Vs.Co$table[,"PValue"], breaks=50, main="Anterior vs Cortex p-value frequency histogram")

# Anterior vs Frontal Cortex
lrt.At.Vs.Fc <- glmLRT(fit, contrast=my.contrasts[,"At.Vs.Fc"])
df_At.Vs.Fc <- summary(decideTests(lrt.At.Vs.Fc))
grid.newpage()
grid.table(df_At.Vs.Fc)
plotMD(lrt.At.Vs.Fc)

# Volcano plot
volcanoData <- cbind(lrt.At.Vs.Fc$table$logFC, -log10(lrt.At.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Frontal_Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.At.Vs.Fc$table[,"PValue"], breaks=50, main="Anterior vs Frontal Cortex p-value frequency histogram")

# Anterior vs Hippocampus
lrt.At.Vs.Hp <- glmLRT(fit, contrast=my.contrasts[,"At.Vs.Hp"])
df_At.Vs.Hp <- summary(decideTests(lrt.At.Vs.Hp))
grid.newpage()
grid.table(df_At.Vs.Hp)
plotMD(lrt.At.Vs.Hp)

# Volcano plot
volcanoData <- cbind(lrt.At.Vs.Hp$table$logFC, -log10(lrt.At.Vs.Hp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.At.Vs.Hp$table[,"PValue"], breaks=50, main="Anterior vs Hippocampus p-value frequency histogram")

# Anterior vs Hypothalamus
lrt.At.Vs.Hy <- glmLRT(fit, contrast=my.contrasts[,"At.Vs.Hy"])
df_At.Vs.Hy <- summary(decideTests(lrt.At.Vs.Hy))
grid.newpage()
grid.table(df_At.Vs.Hy)
plotMD(lrt.At.Vs.Hy)

# Volcano plot
volcanoData <- cbind(lrt.At.Vs.Hy$table$logFC, -log10(lrt.At.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.At.Vs.Hp$table[,"PValue"], breaks=50, main="Anterior vs Hypothalamus p-value frequency histogram")

# Anterior vs Nucleus Accumbens
lrt.At.Vs.Nc <- glmLRT(fit, contrast=my.contrasts[,"At.Vs.Nc"])
df_At.Vs.Nc <- summary(decideTests(lrt.At.Vs.Nc))
grid.newpage()
grid.table(df_At.Vs.Nc)
plotMD(lrt.At.Vs.Nc)

# Volcano plot
volcanoData <- cbind(lrt.At.Vs.Nc$table$logFC, -log10(lrt.At.Vs.Nc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.At.Vs.Nc$table[,"PValue"], breaks=50, main="Anterior vs Nucleus Accumbens p-value frequency histogram")

# Anterior vs Putamen
lrt.At.Vs.Pu <- glmLRT(fit, contrast=my.contrasts[,"At.Vs.Pu"])
df_At.Vs.Pu <- summary(decideTests(lrt.At.Vs.Pu))
grid.newpage()
grid.table(df_At.Vs.Pu)
plotMD(lrt.At.Vs.Pu)

# Volcano plot
volcanoData <- cbind(lrt.At.Vs.Pu$table$logFC, -log10(lrt.At.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.At.Vs.Pu$table[,"PValue"], breaks=50, main="Anterior vs Putamen p-value frequency histogram")

# Anterior vs Spinal Cord
lrt.At.Vs.Sp <- glmLRT(fit, contrast=my.contrasts[,"At.Vs.Sp"])
df_At.Vs.Sp <- summary(decideTests(lrt.At.Vs.Sp))
grid.newpage()
grid.table(df_At.Vs.Sp)
plotMD(lrt.At.Vs.Sp)

# Volcano plot
volcanoData <- cbind(lrt.At.Vs.Sp$table$logFC, -log10(lrt.At.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.At.Vs.Sp$table[,"PValue"], breaks=50, main="Anterior vs Spinal Cord p-value frequency histogram")

# Anterior vs Substantia Nigra
lrt.At.Vs.Sn <- glmLRT(fit, contrast=my.contrasts[,"At.Vs.Sn"])
df_At.Vs.Sn <- summary(decideTests(lrt.At.Vs.Sn))
grid.newpage()
grid.table(df_At.Vs.Sn)
plotMD(lrt.At.Vs.Sn)

# Volcano plot
volcanoData <- cbind(lrt.At.Vs.Sn$table$logFC, -log10(lrt.At.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.At.Vs.Sn$table[,"PValue"], breaks=50, main="Anterior vs Substantia Nigra p-value frequency histogram")

# Caudate vs Cerebellar
lrt.Ca.Vs.Ce <- glmLRT(fit, contrast=my.contrasts[,"Ca.Vs.Ce"])
df_Ca.Vs.Ce <- summary(decideTests(lrt.Ca.Vs.Ce))
grid.newpage()
grid.table(df_Ca.Vs.Ce)
plotMD(lrt.Ca.Vs.Ce)

# Volcano plot
volcanoData <- cbind(lrt.Ca.Vs.Ce$table$logFC, -log10(lrt.Ca.Vs.Ce$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Cerebellar')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ca.Vs.Ce$table[,"PValue"], breaks=50, main="Caudate vs Cerebellar p-value frequency histogram")

# Caudate vs Cerebellum
lrt.Ca.Vs.Cm <- glmLRT(fit, contrast=my.contrasts[,"Ca.Vs.Cm"])
df_Ca.Vs.Cm <- summary(decideTests(lrt.Ca.Vs.Cm))
grid.newpage()
grid.table(df_Ca.Vs.Cm)
plotMD(lrt.At.Vs.Sn)

# Volcano plot
volcanoData <- cbind(lrt.Ca.Vs.Cm$table$logFC, -log10(lrt.Ca.Vs.Cm$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Cerebellum')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ca.Vs.Cm$table[,"PValue"], breaks=50, main="Caudate vs Cerebellum p-value frequency histogram")

# Caudate vs Cortex
lrt.Ca.Vs.Co <- glmLRT(fit, contrast=my.contrasts[,"Ca.Vs.Co"])
df_Ca.Vs.Co <- summary(decideTests(lrt.Ca.Vs.Co))
grid.newpage()
grid.table(df_Ca.Vs.Co)
plotMD(lrt.Ca.Vs.Co)

# Volcano plot
volcanoData <- cbind(lrt.Ca.Vs.Co$table$logFC, -log10(lrt.Ca.Vs.Co$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ca.Vs.Co$table[,"PValue"], breaks=50, main="Caudate vs Cortex p-value frequency histogram")

# Caudate vs Frontal Cortex
lrt.Ca.Vs.Fc <- glmLRT(fit, contrast=my.contrasts[,"Ca.Vs.Fc"])
df_Ca.Vs.Fc <- summary(decideTests(lrt.Ca.Vs.Fc))
grid.newpage()
grid.table(df_Ca.Vs.Fc)
plotMD(lrt.Ca.Vs.Fc)

# Volcano plot
volcanoData <- cbind(lrt.Ca.Vs.Fc$table$logFC, -log10(lrt.Ca.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Frontal_Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ca.Vs.Fc$table[,"PValue"], breaks=50, main="Caudate vs Frontal Cortex p-value frequency histogram")

# Caudate vs Hipocampus
lrt.Ca.Vs.Hp <- glmLRT(fit, contrast=my.contrasts[,"Ca.Vs.Hp"])
df_Ca.Vs.Hp <- summary(decideTests(lrt.Ca.Vs.Hp))
grid.newpage()
grid.table(df_Ca.Vs.Hp)
plotMD(lrt.Ca.Vs.Hp)

# Volcano plot
volcanoData <- cbind(lrt.Ca.Vs.Hp$table$logFC, -log10(lrt.Ca.Vs.Hp[$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ca.Vs.Hp$table[,"PValue"], breaks=50, main="Caudate vs Hippocampus p-value frequency histogram")

# Caudate vs Hypothalamus
lrt.Ca.Vs.Hy <- glmLRT(fit, contrast=my.contrasts[,"Ca.Vs.Hy"])
df_Ca.Vs.Hy <- summary(decideTests(lrt.Ca.Vs.Hy))
grid.newpage()
grid.table(df_Ca.Vs.Hy)
plotMD(lrt.Ca.Vs.Hy)

# Volcano plot
volcanoData <- cbind(lrt.Ca.Vs.Hy$table$logFC, -log10(lrt.Ca.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ca.Vs.Hy$table[,"PValue"], breaks=50, main="Caudate vs Hypothalamus p-value frequency histogram")

# Caudate vs Nucleus Accumbens
lrt.Ca.Vs.Na <- glmLRT(fit, contrast=my.contrasts[,"Ca.Vs.Na"])
df_Ca.Vs.Na <- summary(decideTests(lrt.Ca.Vs.Na))
grid.newpage()
grid.table(df_Ca.Vs.Na)
plotMD(lrt.Ca.Vs.Na)

# Volcano plot
volcanoData <- cbind(lrt.Ca.Vs.Na$table$logFC, -log10(lrt.Ca.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ca.Vs.Na$table[,"PValue"], breaks=50, main="Caudate vs Nucleus Accumbens p-value frequency histogram")

# Caudate vs Putamen
lrt.Ca.Vs.Pu <- glmLRT(fit, contrast=my.contrasts[,"Ca.Vs.Pu"])
df_Ca.Vs.Pu <- summary(decideTests(lrt.Ca.Vs.Pu))
grid.newpage()
grid.table(df_Ca.Vs.Pu)
plotMD(lrt.Ca.Vs.Pu)

# Volcano plot
volcanoData <- cbind(lrt.Ca.Vs.Pu$table$logFC, -log10(lrt.Ca.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ca.Vs.Pu$table[,"PValue"], breaks=50, main="Caudate vs Putamen p-value frequency histogram")

# Caudate vs Spinal Cord
lrt.Ca.Vs.Sp <- glmLRT(fit, contrast=my.contrasts[,"Ca.Vs.Sp"])
df_Ca.Vs.Sp <- summary(decideTests(lrt.Ca.Vs.Sp))
grid.newpage()
grid.table(df_Ca.Vs.Sp)
plotMD(lrt.Ca.Vs.Sp)

# Volcano plot
volcanoData <- cbind(lrt.Ca.Vs.Sp$table$logFC, -log10(lrt.Ca.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ca.Vs.Sp$table[,"PValue"], breaks=50, main="Caudate vs Spinal Cord p-value frequency histogram")

# Caudate vs Substantia Nigra
lrt.Ca.Vs.Sn <- glmLRT(fit, contrast=my.contrasts[,"Ca.Vs.Sn"])
df_Ca.Vs.Sn <- summary(decideTests(lrt.Ca.Vs.Sn))
grid.newpage()
grid.table(df_Ca.Vs.Sn)
plotMD(lrt.Ca.Vs.Sn)

# Volcano plot
volcanoData <- cbind(lrt.Ca.Vs.Sn$table$logFC, -log10(lrt.Ca.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ca.Vs.Sn$table[,"PValue"], breaks=50, main="Caudate vs Substantia Nigra p-value frequency histogram")

# Cerebellar vs Cerebellum
lrt.Ce.Vs.Cm <- glmLRT(fit, contrast=my.contrasts[,"Ce.Vs.Cm"])
df_Ce.Vs.Cm <- summary(decideTests(lrt.Ce.Vs.Cm))
grid.newpage()
grid.table(df_Ce.Vs.Cm)
plotMD(lrt.Ce.Vs.Cm)

# Volcano plot
volcanoData <- cbind(lrt.Ce.Vs.Cm$table$logFC, -log10(lrt.Ce.Vs.Cm$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Cerebellum')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ce.Vs.Cm$table[,"PValue"], breaks=50, main="Cerebellar vs Cerebellum p-value frequency histogram")

# Cerebellar vs Cortex
lrt.Ce.Vs.Co <- glmLRT(fit, contrast=my.contrasts[,"Ce.Vs.Co"])
df_Ce.Vs.Co <- summary(decideTests(lrt.Ce.Vs.Co))
grid.newpage()
grid.table(df_Ce.Vs.Co)
plotMD(lrt.Ce.Vs.Co)

# Volcano plot
volcanoData <- cbind(lrt.Ce.Vs.Co$table$logFC, -log10(lrt.Ce.Vs.Co$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ce.Vs.Co$table[,"PValue"], breaks=50, main="Cerebellar vs Cortex p-value frequency histogram")

# Cerebellar vs Frontal Cortex
lrt.Ce.Vs.Fc <- glmLRT(fit, contrast=my.contrasts[,"Ce.Vs.Fc"])
df_Ce.Vs.Fc <- summary(decideTests(lrt.Ce.Vs.Fc))
grid.newpage()
grid.table(df_Ce.Vs.Fc)
plotMD(lrt.Ce.Vs.Fc)

# Volcano plot
volcanoData <- cbind(lrt.Ce.Vs.Fc$table$logFC, -log10(lrt.Ce.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Frontal_Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ce.Vs.Fc$table[,"PValue"], breaks=50, main="Cerebellar vs Frontal Cortex p-value frequency histogram")

# Cerebellar vs Hippocampus
lrt.Ce.Vs.Hp <- glmLRT(fit, contrast=my.contrasts[,"Ce.Vs.Hp"])
df_Ce.Vs.Hp <- summary(decideTests(lrt.Ce.Vs.Hp))
grid.newpage()
grid.table(df_Ce.Vs.Hp)
plotMD(lrt.Ce.Vs.Hp)

# Volcano plot
volcanoData <- cbind(lrt.Ce.Vs.Hp$table$logFC, -log10(lrt.Ce.Vs.Hp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ce.Vs.Hp$table[,"PValue"], breaks=50, main="Cerebellar vs Hippocampus p-value frequency histogram")

# Cerebellar vs Hypothalamus
lrt.Ce.Vs.Hy <- glmLRT(fit, contrast=my.contrasts[,"Ce.Vs.Hy"])
df_Ce.Vs.Hy <- summary(decideTests(lrt.Ce.Vs.Hy))
grid.newpage()
grid.table(df_Ce.Vs.Hy)
plotMD(lrt.Ce.Vs.Hy)

# Volcano plot
volcanoData <- cbind(lrt.Ce.Vs.Hy$table$logFC, -log10(lrt.Ce.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ce.Vs.Hy$table[,"PValue"], breaks=50, main="Cerebellar vs Hypothalamus p-value frequency histogram")

# Cerebellar vs Nucleus Accumbens
lrt.Ce.Vs.Na <- glmLRT(fit, contrast=my.contrasts[,"Ce.Vs.Na"])
df_Ce.Vs.Na <- summary(decideTests(lrt.Ce.Vs.Na))
grid.newpage()
grid.table(df_Ce.Vs.Na)
plotMD(lrt.Ce.Vs.Na)

# Volcano plot
volcanoData <- cbind(lrt.Ce.Vs.Na$table$logFC, -log10(lrt.Ce.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ce.Vs.Na$table[,"PValue"], breaks=50, main="Cerebellar vs Nucleus Accumbens p-value frequency histogram")

# Cerebellar vs Putamen
lrt.Ce.Vs.Pu <- glmLRT(fit, contrast=my.contrasts[,"Ce.Vs.Pu"])
df_Ce.Vs.Pu <- summary(decideTests(lrt.Ce.Vs.Pu))
grid.newpage()
grid.table(df_Ce.Vs.Pu)
plotMD(lrt.Ce.Vs.Pu)

# Volcano plot
volcanoData <- cbind(lrt.Ce.Vs.Pu$table$logFC, -log10(lrt.Ce.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ce.Vs.Pu$table[,"PValue"], breaks=50, main="Cerebellar vs Putamen p-value frequency histogram")

# Cerebellar vs Spinal Cord
lrt.Ce.Vs.Sp <- glmLRT(fit, contrast=my.contrasts[,"Ce.Vs.Sp"])
df_Ce.Vs.Sp <- summary(decideTests(lrt.Ce.Vs.Sp))
grid.newpage()
grid.table(df_Ce.Vs.Sp)
plotMD(lrt.Ce.Vs.Sp)

# Volcano plot
volcanoData <- cbind(lrt.Ce.Vs.Sp$table$logFC, -log10(lrt.Ce.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ce.Vs.Sp$table[,"PValue"], breaks=50, main="Cerebellar vs Spinal Cord p-value frequency histogram")

# Cerebellar vs Substantia Nigra
lrt.Ce.Vs.Sn <- glmLRT(fit, contrast=my.contrasts[,"Ce.Vs.Sn"])
df_Ce.Vs.Sn <- summary(decideTests(lrt.Ce.Vs.Sn))
grid.newpage()
grid.table(df_Ce.Vs.Sn)
plotMD(lrt.Ce.Vs.Sn)

# Volcano plot
volcanoData <- cbind(lrt.Ce.Vs.Sn$table$logFC, -log10(lrt.Ce.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Ce.Vs.Sn$table[,"PValue"], breaks=50, main="Cerebellar vs Substantia Nigra p-value frequency histogram")

# Cortex vs Frontal Cortex
lrt.Co.Vs.Fc <- glmLRT(fit, contrast=my.contrasts[,"Co.Vs.Fc"])
df_Co.Vs.Fc <- summary(decideTests(lrt.Co.Vs.Fc))
grid.newpage()
grid.table(df_Co.Vs.Fc)
plotMD(lrt.Co.Vs.Fc)

# Volcano plot
volcanoData <- cbind(lrt.Co.Vs.Fc$table$logFC, -log10(lrt.Co.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortexr-1*Frontal_Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Co.Vs.Fc$table[,"PValue"], breaks=50, main="Cortex vs Frontal Cortex p-value frequency histogram")

# Cortex vs Hippocampus
lrt.Co.Vs.Hp <- glmLRT(fit, contrast=my.contrasts[,"Co.Vs.Hp"])
df_Co.Vs.Hp <- summary(decideTests(lrt.Co.Vs.Hp))
grid.newpage()
grid.table(df_Co.Vs.Hp)
plotMD(lrt.Co.Vs.Hp)

# Volcano plot
volcanoData <- cbind(lrt.Co.Vs.Hp$table$logFC, -log10(lrt.Co.Vs.Hp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Co.Vs.Hp$table[,"PValue"], breaks=50, main="Cortex vs Hippocampus p-value frequency histogram")

# Cortex vs Hypothalamus
lrt.Co.Vs.Hy <- glmLRT(fit, contrast=my.contrasts[,"Co.Vs.Hy"])
df_Co.Vs.Hy <- summary(decideTests(lrt.Co.Vs.Hy))
grid.newpage()
grid.table(df_Co.Vs.Hy)
plotMD(lrt.Co.Vs.Hy)

# Volcano plot
volcanoData <- cbind(lrt.Co.Vs.Fc$table$logFC, -log10(lrt.Co.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Co.Vs.Hy$table[,"PValue"], breaks=50, main="Cortex vs Hypothalamus p-value frequency histogram")

# Cortex vs Nucleus Accumbens
lrt.Co.Vs.Na <- glmLRT(fit, contrast=my.contrasts[,"Co.Vs.Na"])
df_Co.Vs.Na <- summary(decideTests(lrt.Co.Vs.Na))
grid.newpage()
grid.table(df_Co.Vs.Na)
plotMD(lrt.Co.Vs.Na)

# Volcano plot
volcanoData <- cbind(lrt.Co.Vs.Na$table$logFC, -log10(lrt.Co.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Co.Vs.Na$table[,"PValue"], breaks=50, main="Cortex vs Nucleus Accumbens p-value frequency histogram")

# Cortex vs Putamen
lrt.Co.Vs.Pu <- glmLRT(fit, contrast=my.contrasts[,"Co.Vs.Pu"])
df_Co.Vs.Pu <- summary(decideTests(lrt.Co.Vs.Pu))
grid.newpage()
grid.table(df_Co.Vs.Pu)
plotMD(lrt.Co.Vs.Pu)

# Volcano plot
volcanoData <- cbind(lrt.Co.Vs.Pu$table$logFC, -log10(lrt.Co.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Co.Vs.Pu$table[,"PValue"], breaks=50, main="Cortex vs Putamen p-value frequency histogram")

# Cortex vs Spinal Cord
lrt.Co.Vs.Sp <- glmLRT(fit, contrast=my.contrasts[,"Co.Vs.Sp"])
df_Co.Vs.Sp <- summary(decideTests(lrt.Co.Vs.Sp))
grid.newpage()
grid.table(df_Co.Vs.Sp)
plotMD(lrt.Co.Vs.Sp)

# Volcano plot
volcanoData <- cbind(lrt.Co.Vs.Sp$table$logFC, -log10(lrt.Co.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Co.Vs.Sp$table[,"PValue"], breaks=50, main="Cortex vs Spinal Cord p-value frequency histogram")

# Cortex vs Substantia Nigra
lrt.Co.Vs.Sn <- glmLRT(fit, contrast=my.contrasts[,"Co.Vs.Sn"])
df_Co.Vs.Sn <- summary(decideTests(lrt.Co.Vs.Sn))
grid.newpage()
grid.table(df_Co.Vs.Sn)
plotMD(lrt.Co.Vs.Sn)

# Volcano plot
volcanoData <- cbind(lrt.Co.Vs.Sn$table$logFC, -log10(lrt.Co.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Co.Vs.Sn$table[,"PValue"], breaks=50, main="Cortex vs Substantia Nigra p-value frequency histogram")

# Frontal Cortex vs Hippocampus
lrt.Fc.Vs.Hp <- glmLRT(fit, contrast=my.contrasts[,"Fc.Vs.Hp"])
df_Fc.Vs.Hp <- summary(decideTests(lrt.Fc.Vs.Hp))
grid.newpage()
grid.table(df_Fc.Vs.Hp)
plotMD(lrt.Fc.Vs.Hp)

# Volcano plot
volcanoData <- cbind(lrt.Fc.Vs.Hp$table$logFC, -log10(lrt.Fc.Vs.Hp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Fc.Vs.Hp$table[,"PValue"], breaks=50, main="Frontal Cortex vs Hippocampus p-value frequency histogram")

# Frontal Cortex vs Hypothalamus 
lrt.Fc.Vs.Hy <- glmLRT(fit, contrast=my.contrasts[,"Fc.Vs.Hy"])
df_Fc.Vs.Hy <- summary(decideTests(lrt.Fc.Vs.Hy))
grid.newpage()
grid.table(df_Fc.Vs.Hy)
plotMD(lrt.Fc.Vs.Hy)

# Volcano plot
volcanoData <- cbind(lrt.Fc.Vs.Hy$table$logFC, -log10(lrt.Fc.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Fc.Vs.Hy$table[,"PValue"], breaks=50, main="Frontal Cortex vs Hypothalamus p-value frequency histogram")

# Frontal Cortex vs Nucleus Accumbens
lrt.Fc.Vs.Na <- glmLRT(fit, contrast=my.contrasts[,"Fc.Vs.Na"])
df_Fc.Vs.Na <- summary(decideTests(lrt.Fc.Vs.Na))
grid.newpage()
grid.table(df_Fc.Vs.Na)
plotMD(lrt.Fc.Vs.Na)

# Volcano plot
volcanoData <- cbind(lrt.Fc.Vs.Na$table$logFC, -log10(lrt.Fc.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Fc.Vs.Na$table[,"PValue"], breaks=50, main="Frontal Cortex vs Nucleus Accumbens p-value frequency histogram")

# Frontal Cortex vs Putamen
lrt.Fc.Vs.Pu <- glmLRT(fit, contrast=my.contrasts[,"Fc.Vs.Pu"])
df_Fc.Vs.Pu <- summary(decideTests(lrt.Fc.Vs.Pu))
grid.newpage()
grid.table(df_Fc.Vs.Pu)
plotMD(lrt.Fc.Vs.Pu)

# Volcano plot
volcanoData <- cbind(lrt.Fc.Vs.Pu$table$logFC, -log10(lrt.Fc.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Fc.Vs.Pu$table[,"PValue"], breaks=50, main="Frontal Cortex vs Putamen p-value frequency histogram")

# Frontal Corte vs Spinal Cord
lrt.Fc.Vs.Sp <- glmLRT(fit, contrast=my.contrasts[,"Fc.Vs.Sp"])
df_Fc.Vs.Sp <- summary(decideTests(lrt.Fc.Vs.Sp))
grid.newpage()
grid.table(df_Fc.Vs.Sp)
plotMD(lrt.Fc.Vs.Sp)

# Volcano plot
volcanoData <- cbind(lrt.Fc.Vs.Sp$table$logFC, -log10(lrt.Fc.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Fc.Vs.Sp$table[,"PValue"], breaks=50, main="Frontal Cortex vs Spinal Cord p-value frequency histogram")

# Frontal Cortex vs Substantia Nigra
lrt.Fc.Vs.Sn <- glmLRT(fit, contrast=my.contrasts[,"Fc.Vs.Sn"])
df_Fc.Vs.Sn <- summary(decideTests(lrt.Fc.Vs.Sn))
grid.newpage()
grid.table(df_Fc.Vs.Sn)
plotMD(lrt.Fc.Vs.Sn)

# Volcano plot
volcanoData <- cbind(lrt.Fc.Vs.Sn$table$logFC, -log10(lrt.Fc.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Fc.Vs.Sn$table[,"PValue"], breaks=50, main="Frontal Cortex vs Substantia Nigra p-value frequency histogram")

# Hippocampus vs Hypothalamus
lrt.Hp.Vs.Hy <- glmLRT(fit, contrast=my.contrasts[,"Hp.Vs.Hy"])
df_Hp.Vs.Hy <- summary(decideTests(lrt.Hp.Vs.Hy))
grid.newpage()
grid.table(df_Hp.Vs.Hy)
plotMD(lrt.Hp.Vs.Hy)

# Volcano plot
volcanoData <- cbind(lrt.Hp.Vs.Hy$table$logFC, -log10(lrt.Hp.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Hp.Vs.Hy$table[,"PValue"], breaks=50, main="Hippocampus vs Hypothalamus p-value frequency histogram")

# Hippocampus vs Nucleus Accumbens
lrt.Hp.Vs.Na <- glmLRT(fit, contrast=my.contrasts[,"Hp.Vs.Na"])
df_Hp.Vs.Na <- summary(decideTests(lrt.Hp.Vs.Na))
grid.newpage()
grid.table(df_Hp.Vs.Na)
plotMD(lrt.Hp.Vs.Na)

# Volcano plot
volcanoData <- cbind(lrt.Hp.Vs.Na$table$logFC, -log10(lrt.Hp.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Hp.Vs.Na$table[,"PValue"], breaks=50, main="Hippocampus vs Nucleus Accumbens p-value frequency histogram")

# Hippocampus vs Putamen
lrt.Hp.Vs.Pu <- glmLRT(fit, contrast=my.contrasts[,"Hp.Vs.Pu"])
df_Hp.Vs.Pu <- summary(decideTests(lrt.Hp.Vs.Pu))
grid.newpage()
grid.table(df_Hp.Vs.Pu)
plotMD(lrt.Hp.Vs.Pu)

# Volcano plot
volcanoData <- cbind(lrt.Hp.Vs.Pu$table$logFC, -log10(lrt.Hp.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Hp.Vs.Pu$table[,"PValue"], breaks=50, main="Hippocampus vs Putamen p-value frequency histogram")

# Hippocampus vs Spinal Cord
lrt.Hp.Vs.Sp <- glmLRT(fit, contrast=my.contrasts[,"Hp.Vs.Sp"])
df_Hp.Vs.Sp <- summary(decideTests(lrt.Hp.Vs.Sp))
grid.newpage()
grid.table(df_Hp.Vs.Sp)
plotMD(lrt.Hp.Vs.Sp)

# Volcano plot
volcanoData <- cbind(lrt.Hp.Vs.Sp$table$logFC, -log10(lrt.Hp.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Hp.Vs.Sp$table[,"PValue"], breaks=50, main="Hippocampus vs Spinal Cord p-value frequency histogram")

# Hippocampus vs Substantia Nigra 
lrt.Fc.Vs.Sn <- glmLRT(fit, contrast=my.contrasts[,"Fc.Vs.Sn"])
df_Hp.Vs.Sn <- summary(decideTests(lrt.Hp.Vs.Sn))
grid.newpage()
grid.table(df_Hp.Vs.Sn)
plotMD(lrt.Hp.Vs.Sn)

# Volcano plot
volcanoData <- cbind(lrt.Hp.Vs.Sn$table$logFC, -log10(lrt.Hp.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Hp.Vs.Sn$table[,"PValue"], breaks=50, main="Hippocampus vs Substantia Nigra p-value frequency histogram")

# Hypothalamus vs Nucleus Accumbens
lrt.Hy.Vs.Na <- glmLRT(fit, contrast=my.contrasts[,"Hy.Vs.Na"])
df_Hy.Vs.Na <- summary(decideTests(lrt.Hy.Vs.Na))
grid.newpage()
grid.table(df_Hy.Vs.Na)
plotMD(lrt.Hy.Vs.Na)

# Volcano plot
volcanoData <- cbind(lrt.Hy.Vs.Na$table$logFC, -log10(lrt.Hy.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hypothalamus-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Hy.Vs.Na$table[,"PValue"], breaks=50, main="Hypothalamus vs Nucleus Accumbens p-value frequency histogram")

# Hypothalamus vs Putamen
lrt.Hy.Vs.Pu <- glmLRT(fit, contrast=my.contrasts[,"Hy.Vs.Pu"])
df_Hy.Vs.Pu <- summary(decideTests(lrt.Hy.Vs.Pu))
grid.newpage()
grid.table(df_Hy.Vs.Pu)
plotMD(lrt.Hy.Vs.Pu)

# Volcano plot
volcanoData <- cbind(lrt.Hy.Vs.Pu$table$logFC, -log10(lrt.Hy.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hypothalamus-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Hy.Vs.Pu$table[,"PValue"], breaks=50, main="Hypothalamus vs Putamen p-value frequency histogram")

# Hypothalamus vs Spinal Cord
lrt.Hy.Vs.Sp <- glmLRT(fit, contrast=my.contrasts[,"Hy.Vs.Sp"])
df_Hy.Vs.Sp <- summary(decideTests(lrt.Hy.Vs.Sp))
grid.newpage()
grid.table(df_Hy.Vs.Sp)
plotMD(lrt.Hy.Vs.Sp)

# Volcano plot
volcanoData <- cbind(lrt.Hy.Vs.Sp$table$logFC, -log10(lrt.Hy.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hypothalamus-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Hy.Vs.Sp$table[,"PValue"], breaks=50, main="Hypothalamus vs Spinal Cord p-value frequency histogram")

# Hypothalamus Substantia Nigra
lrt.Hy.Vs.Sn <- glmLRT(fit, contrast=my.contrasts[,"Hy.Vs.Sn"])
df_Hy.Vs.Sn <- summary(decideTests(lrt.Hy.Vs.Sn))
grid.newpage()
grid.table(df_Hy.Vs.Sn)
plotMD(lrt.Hy.Vs.Sn)

# Volcano plot
volcanoData <- cbind(lrt.Hy.Vs.Sn$table$logFC, -log10(lrt.Hy.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hypothalamus-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Hy.Vs.Sn$table[,"PValue"], breaks=50, main="Hypothalamus vs Substantia Nigra p-value frequency histogram")

# Nucleus Accumbens vs Putamen
lrt.Na.Vs.Pu <- glmLRT(fit, contrast=my.contrasts[,"Na.Vs.Pu"])
df_Na.Vs.Pu <- summary(decideTests(lrt.Na.Vs.Pu))
grid.newpage()
grid.table(df_Na.Vs.Pu)
plotMD(lrt.Na.Vs.Pu)

# Volcano plot
volcanoData <- cbind(lrt.Na.Vs.Pu$table$logFC, -log10(lrt.Na.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Nucleus_Accumbens-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Na.Vs.Pu$table[,"PValue"], breaks=50, main="Nucleus Accumbens vs Putamen p-value frequency histogram")

# Nucleus Accumbens vs Spinal Cord
lrt.Na.Vs.Sp <- glmLRT(fit, contrast=my.contrasts[,"Na.Vs.Sp"])
df_Na.Vs.Sp <- summary(decideTests(lrt.Na.Vs.Sp))
grid.newpage()
grid.table(df_Na.Vs.Sp)
plotMD(lrt.Na.Vs.Sp)

# Volcano plot
volcanoData <- cbind(lrt.Na.Vs.Sp$table$logFC, -log10(lrt.Na.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Nucleus_Accumbens-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Na.Vs.Sp$table[,"PValue"], breaks=50, main="Nucleus Accumbens vs Spinal Cord p-value frequency histogram")

# Nucleus Accumbens vs Substantia Nigra
lrt.Na.Vs.Sn <- glmLRT(fit, contrast=my.contrasts[,"Na.Vs.Sn"])
df_Na.Vs.Sn <- summary(decideTests(lrt.Na.Vs.Sn))
grid.newpage()
grid.table(df_na.Vs.Sn)
plotMD(lrt.Na.Vs.Sn)

# Volcano plot
volcanoData <- cbind(lrt.Na.Vs.Sn$table$logFC, -log10(lrt.Na.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Nucleus_Accumbens-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Na.Vs.Sn$table[,"PValue"], breaks=50, main="Nucleus Accumbens vs Substantia Nigra p-value frequency histogram")

# Putamen vs Spinal Cord
lrt.Pu.Vs.Sp <- glmLRT(fit, contrast=my.contrasts[,"Pu.Vs.Sp"])
df_Pu.Vs.Sp <- summary(decideTests(lrt.Pu.Vs.Sp))
grid.newpage()
grid.table(df_Pu.Vs.Sp)
plotMD(lrt.Pu.Vs.Sp)

# Volcano plot
volcanoData <- cbind(lrt.Pu.Vs.Sp$table$logFC, -log10(lrt.Pu.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Putamen-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Pu.Vs.Sp$table[,"PValue"], breaks=50, main="Putamen vs Spinal Cord p-value frequency histogram")

# Putamen vs Substantia Nigra
lrt.Pu.Vs.Sn <- glmLRT(fit, contrast=my.contrasts[,"Pu.Vs.Sn"])
df_Pu.Vs.Sn <- summary(decideTests(lrt.Pu.Vs.Sn))
grid.newpage()
grid.table(df_Pu.Vs.Sn)
plotMD(lrt.Pu.Vs.Sn)

# Volcano plot
volcanoData <- cbind(lrt.Pu.Vs.Sn$table$logFC, -log10(lrt.Pu.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Putamen-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Pu.Vs.Sn$table[,"PValue"], breaks=50, main="Putamen vs Substantia Nigra p-value frequency histogram")

# Spinal Cord vs Substantia Nigra
lrt.Sp.Vs.Sn <- glmLRT(fit, contrast=my.contrasts[,"Sp.Vs.Sn"])
df_Sp.Vs.Sn <- summary(decideTests(lrt.Sp.Vs.Sn))
grid.newpage()
grid.table(df_Sp.Vs.Sn)
plotMD(lrt.Sp.Vs.Sn)

# Volcano plot
volcanoData <- cbind(lrt.Sp.Vs.Sn$table$logFC, -log10(lrt.Sp.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Spinal_Cord-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Sp.Vs.Sn$table[,"PValue"], breaks=50, main="Spinal Cord vs Substantia Nigra p-value frequency histogram")

dev.off()

