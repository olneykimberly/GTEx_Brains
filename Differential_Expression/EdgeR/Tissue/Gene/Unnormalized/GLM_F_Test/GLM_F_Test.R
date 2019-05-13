# This script looks at differential gene expression between each brain tissue type.

METADATA <- '/scratch/mjpete11/GTEx/Toy_Metadata.csv'
COUNT_MATRIX <- '/scratch/mjpete11/GTEx/Data_Exploration/Toy_Count_Matrix.tsv'
PLOT_DIR =  '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Tissue/Gene/Normalized/GLM_F_Test/Plots/'
FILE_NAME = 'Unnormalized_GLM_F_Test.pdf'

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

# Test for DEX genes using quasi-liklihood method.
fit <- glmQLFit(y, design, robust=TRUE)

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
                              levels = design)


# Plot
setwd(PLOT_DIR)

pdf(FILE_NAME)

# Amygdala vs Anterior
qlf.Am.Vs.At <- glmQLFTest(fit, contrast=my.contrasts[,"Am.Vs.At"])
df_Am.Vs.At <- summary(decideTests(qlf.Am.Vs.At))
grid.newpage()
grid.table(df_Am.Vs.At)
plotMD(qlf.Am.Vs.At)

# Volcano plot
volcanoData <- cbind(qlf.Am.Vs.At$table$logFC, -log10(qlf.Am.Vs.At$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Anterior')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Am.Vs.At$table[,"PValue"], breaks=50, main="Amygdala vs Anterior p-value frequency histogram")

# Amygdala vs Caduate
qlf.Am.Vs.Ca <- glmQLFTest(fit, contrast=my.contrasts[,"Am.Vs.Ca"])
df_Am.Vs.Ca <- summary(decideTests(qlf.Am.Vs.Ca))
grid.newpage()
grid.table(df_Am.Vs.Ca)
plotMD(qlf.Am.Vs.Ca)

# Volcano plot
volcanoData <- cbind(qlf.Am.Vs.Ca$table$logFC, -log10(qlf.Am.Vs.Ca$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Caudate')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Am.Vs.Ca$table[,"PValue"], breaks=50, main="Amygdala vs Caudate p-value frequency histogram")

# Amygdala vs Cerebellar
qlf.Am.Vs.Ce <- glmQLFTest(fit, contrast=my.contrasts[,"Am.Vs.Ce"])
df_Am.Vs.Ce <- summary(decideTests(qlf.Am.Vs.Ce))
grid.newpage()
grid.table(df_Am.Vs.Ce)
plotMD(qlf.Am.Vs.Ce)

# Volcano plot
volcanoData <- cbind(qlf.Am.Vs.Ce$table$logFC, -log10(qlf.Am.Vs.Ce$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Cerebellar')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Am.Vs.Ce$table[,"PValue"], breaks=50, main="Amygdala vs Cerebellar p-value frequency histogram")

# Amygdala vs Cerebellum
qlf.Am.Vs.Cm <- glmQLFTest(fit, contrast=my.contrasts[,"Am.Vs.Cm"])
df_Am.Vs.Cm <- summary(decideTests(qlf.Am.Vs.Cm))
grid.newpage()
grid.table(df_Am.Vs.Cm)
plotMD(qlf.Am.Vs.Cm)

# Volcano plot
volcanoData <- cbind(qlf.Am.Vs.Cm$table$logFC, -log10(qlf.Am.Vs.Cm$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Cerebellum')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Am.Vs.Cm$table[,"PValue"], breaks=50, main="Amygdala vs Cerebellum p-value frequency histogram")

# Amygdala vs Cortex
qlf.Am.Vs.Co <- glmQLFTest(fit, contrast=my.contrasts[,"Am.Vs.Co"])
df_Am.Vs.Co <- summary(decideTests(qlf.Am.Vs.Co))
grid.newpage()
grid.table(df_Am.Vs.Co)
plotMD(qlf.Am.Vs.Co)

# Volcano plot
volcanoData <- cbind(qlf.Am.Vs.Co$table$logFC, -log10(qlf.Am.Vs.Co$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Am.Vs.Co$table[,"PValue"], breaks=50, main="Amygdala vs Cortex p-value frequency histogram")

# Amygdala vs Frontal Cortex
qlf.Am.Vs.Fc <- glmQLFTest(fit, contrast=my.contrasts[,"Am.Vs.Fc"])
df_Am.Vs.Fc <- summary(decideTests(qlf.Am.Vs.Fc))
grid.newpage()
grid.table(df_Am.Vs.Fc)
plotMD(qlf.Am.Vs.Fc)

# Volcano plot
volcanoData <- cbind(qlf.Am.Vs.Fc$table$logFC, -log10(qlf.Am.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Frontal_Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Am.Vs.Fc$table[,"PValue"], breaks=50, main="Amygdala vs Frontal Cortex p-value frequency histogram")

# Amygdala vs Hippocampus
qlf.Am.Vs.Hp <- glmQLFTest(fit, contrast=my.contrasts[,"Am.Vs.Hp"])
summary(decideTests(qlf.Am.Vs.Hp))
plotMD(qlf.Am.Vs.Hp)

# Volcano plot
volcanoData <- cbind(qlf.Am.Vs.Hp$table$logFC, -log10(qlf.Am.Vs.Hp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Am.Vs.Hp$table[,"PValue"], breaks=50, main="Amygdala vs Hippocampus p-value frequency histogram")

# Amygdala vs Hypothalamus
qlf.Am.Vs.Hy <- glmQLFTest(fit, contrast=my.contrasts[,"Am.Vs.Hy"])
df_Am.Vs.Hy <- summary(decideTests(qlf.Am.Vs.Hy))
grid.newpage()
grid.table(df_Am.Vs.Hy)
plotMD(qlf.Am.Vs.Hy)

# Volcano plot
volcanoData <- cbind(qlf.Am.Vs.Hy$table$logFC, -log10(qlf.Am.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Am.Vs.Hy$table[,"PValue"], breaks=50, main="Amygdala vs Hypothalamus p-value frequency histogram")

# Amygdala vs Nucleus Accumbens
qlf.Am.Vs.Nc <- glmQLFTest(fit, contrast=my.contrasts[,"Am.Vs.Nc"])
df_Am.Vs.Nc <- summary(decideTests(qlf.Am.Vs.Nc))
grid.newpage()
grid.table(df_Am.Vs.Nc)
plotMD(qlf.Am.Vs.Nc)

# Volcano plot
volcanoData <- cbind(qlf.Am.Vs.Nc$table$logFC, -log10(qlf.Am.Vs.Nc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Am.Vs.Nc$table[,"PValue"], breaks=50, main="Amygdala vs Nucleus Accumbens p-value frequency histogram")

# Amygdala vs Putamen
qlf.Am.Vs.Pu <- glmQLFTest(fit, contrast=my.contrasts[,"Am.Vs.Pu"])
df_Am.Vs.Pu <- summary(decideTests(qlf.Am.Vs.Pu))
grid.newpage()
grid.table(df_Am.Vs.Pu)
plotMD(qlf.Am.Vs.Pu)

# Volcano plot
volcanoData <- cbind(qlf.Am.Vs.Pu$table$logFC, -log10(qlf.Am.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Am.Vs.Pu$table[,"PValue"], breaks=50, main="Amygdala vs Putamen p-value frequency histogram")

# Amygdala vs Spinal Cord
qlf.Am.Vs.Sp <- glmQLFTest(fit, contrast=my.contrasts[,"Am.Vs.Sp"])
df_Am.Vs.Sp <- summary(decideTests(qlf.Am.Vs.Sp))
grid.newpage()
grid.table(df_Am.Vs.Sp)
plotMD(qlf.Am.Vs.Sp)

# Volcano plot
volcanoData <- cbind(qlf.Am.Vs.Sp$table$logFC, -log10(qlf.Am.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Am.Vs.Sp$table[,"PValue"], breaks=50, main="Amygdala vs Spinal Cord p-value frequency histogram")

# Amygdala vs Substantia Nigra
qlf.Am.Vs.Sn <- glmQLFTest(fit, contrast=my.contrasts[,"Am.Vs.Sn"])
df_Am.Vs.Sn <- summary(decideTests(qlf.Am.Vs.Sn))
grid.newpage()
grid.table(df_Am.Vs.Sn)
plotMD(qlf.Am.Vs.Sn)

# Volcano plot
volcanoData <- cbind(qlf.Am.Vs.Sn$table$logFC, -log10(qlf.Am.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Amygdala-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Am.Vs.Sn$table[,"PValue"], breaks=50, main="Amygdala vs Substantia_Nigra p-value frequency histogram")

# Anterior vs Caudate
qlf.At.Vs.Ca <- glmQLFTest(fit, contrast=my.contrasts[,"At.Vs.Ca"])
df_At.Vs.Ca <- summary(decideTests(qlf.At.Vs.Ca))
grid.newpage()
grid.table(df_At.Vs.Ca)
plotMD(qlf.At.Vs.Ca)

# Volcano plot
volcanoData <- cbind(qlf.At.Vs.Ca$table$logFC, -log10(qlf.At.Vs.Ca$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Caudate')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.At.Vs.Ca$table[,"PValue"], breaks=50, main="Anterior vs Caudate p-value frequency histogram")

# Anterior vs Cerebellar
qlf.At.Vs.Ce <- glmQLFTest(fit, contrast=my.contrasts[,"At.Vs.Ce"])
df_At.Vs.Ce <- summary(decideTests(qlf.At.Vs.Ce))
grid.newpage()
grid.table(df_At.Vs.Ce)
plotMD(qlf.At.Vs.Ce)

# Volcano plot
volcanoData <- cbind(qlf.At.Vs.Ce$table$logFC, -log10(qlf.At.Vs.Ce$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Cerebellar')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.At.Vs.Ce$table[,"PValue"], breaks=50, main="Anterior vs Cerebellar p-value frequency histogram")

# Anterior vs Cerebellum
qlf.At.Vs.Cm <- glmQLFTest(fit, contrast=my.contrasts[,"At.Vs.Cm"])
df_At.Vs.Cm <- summary(decideTests(qlf.At.Vs.Cm))
grid.newpage()
grid.table(df_At.Vs.Cm)
plotMD(qlf.At.Vs.Cm)

# Volcano plot
volcanoData <- cbind(qlf.At.Vs.Cm$table$logFC, -log10(qlf.At.Vs.Cm$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Cerebellum')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.At.Vs.Cm$table[,"PValue"], breaks=50, main="Anterior vs Cerebellum p-value frequency histogram")

# Anterior vs Cortex
qlf.At.Vs.Co <- glmQLFTest(fit, contrast=my.contrasts[,"At.Vs.Co"])
df_At.Vs.Co <- summary(decideTests(qlf.At.Vs.Co))
grid.newpage()
grid.table(df_At.Vs.Co)
plotMD(qlf.At.Vs.Co)

# Volcano plot
volcanoData <- cbind(qlf.At.Vs.Co$table$logFC, -log10(qlf.At.Vs.Co$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.At.Vs.Co$table[,"PValue"], breaks=50, main="Anterior vs Cortex p-value frequency histogram")

# Anterior vs Frontal Cortex
qlf.At.Vs.Fc <- glmQLFTest(fit, contrast=my.contrasts[,"At.Vs.Fc"])
df_At.Vs.Fc <- summary(decideTests(qlf.At.Vs.Fc))
grid.newpage()
grid.table(df_At.Vs.Fc)
plotMD(qlf.At.Vs.Fc)

# Volcano plot
volcanoData <- cbind(qlf.At.Vs.Fc$table$logFC, -log10(qlf.At.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Frontal_Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.At.Vs.Fc$table[,"PValue"], breaks=50, main="Anterior vs Frontal Cortex p-value frequency histogram")

# Anterior vs Hippocampus
qlf.At.Vs.Hp <- glmQLFTest(fit, contrast=my.contrasts[,"At.Vs.Hp"])
df_At.Vs.Hp <- summary(decideTests(qlf.At.Vs.Hp))
grid.newpage()
grid.table(df_At.Vs.Hp)
plotMD(qlf.At.Vs.Hp)

# Volcano plot
volcanoData <- cbind(qlf.At.Vs.Hp$table$logFC, -log10(qlf.At.Vs.Hp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.At.Vs.Hp$table[,"PValue"], breaks=50, main="Anterior vs Hippocampus p-value frequency histogram")

# Anterior vs Hypothalamus
qlf.At.Vs.Hy <- glmQLFTest(fit, contrast=my.contrasts[,"At.Vs.Hy"])
df_At.Vs.Hy <- summary(decideTests(qlf.At.Vs.Hy))
grid.newpage()
grid.table(df_At.Vs.Hy)
plotMD(qlf.At.Vs.Hy)

# Volcano plot
volcanoData <- cbind(qlf.At.Vs.Hy$table$logFC, -log10(qlf.At.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.At.Vs.Hp$table[,"PValue"], breaks=50, main="Anterior vs Hypothalamus p-value frequency histogram")

# Anterior vs Nucleus Accumbens
qlf.At.Vs.Nc <- glmQLFTest(fit, contrast=my.contrasts[,"At.Vs.Nc"])
df_At.Vs.Nc <- summary(decideTests(qlf.At.Vs.Nc))
grid.newpage()
grid.table(df_At.Vs.Nc)
plotMD(qlf.At.Vs.Nc)

# Volcano plot
volcanoData <- cbind(qlf.At.Vs.Nc$table$logFC, -log10(qlf.At.Vs.Nc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.At.Vs.Nc$table[,"PValue"], breaks=50, main="Anterior vs Nucleus Accumbens p-value frequency histogram")

# Anterior vs Putamen
qlf.At.Vs.Pu <- glmQLFTest(fit, contrast=my.contrasts[,"At.Vs.Pu"])
df_At.Vs.Pu <- summary(decideTests(qlf.At.Vs.Pu))
grid.newpage()
grid.table(df_At.Vs.Pu)
plotMD(qlf.At.Vs.Pu)

# Volcano plot
volcanoData <- cbind(qlf.At.Vs.Pu$table$logFC, -log10(qlf.At.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.At.Vs.Pu$table[,"PValue"], breaks=50, main="Anterior vs Putamen p-value frequency histogram")

# Anterior vs Spinal Cord
qlf.At.Vs.Sp <- glmQLFTest(fit, contrast=my.contrasts[,"At.Vs.Sp"])
df_At.Vs.Sp <- summary(decideTests(qlf.At.Vs.Sp))
grid.newpage()
grid.table(df_At.Vs.Sp)
plotMD(qlf.At.Vs.Sp)

# Volcano plot
volcanoData <- cbind(qlf.At.Vs.Sp$table$logFC, -log10(qlf.At.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.At.Vs.Sp$table[,"PValue"], breaks=50, main="Anterior vs Spinal Cord p-value frequency histogram")

# Anterior vs Substantia Nigra
qlf.At.Vs.Sn <- glmQLFTest(fit, contrast=my.contrasts[,"At.Vs.Sn"])
df_At.Vs.Sn <- summary(decideTests(qlf.At.Vs.Sn))
grid.newpage()
grid.table(df_At.Vs.Sn)
plotMD(qlf.At.Vs.Sn)

# Volcano plot
volcanoData <- cbind(qlf.At.Vs.Sn$table$logFC, -log10(qlf.At.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Anterior-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.At.Vs.Sn$table[,"PValue"], breaks=50, main="Anterior vs Substantia Nigra p-value frequency histogram")

# Caudate vs Cerebellar
qlf.Ca.Vs.Ce <- glmQLFTest(fit, contrast=my.contrasts[,"Ca.Vs.Ce"])
df_Ca.Vs.Ce <- summary(decideTests(qlf.Ca.Vs.Ce))
grid.newpage()
grid.table(df_Ca.Vs.Ce)
plotMD(qlf.Ca.Vs.Ce)

# Volcano plot
volcanoData <- cbind(qlf.Ca.Vs.Ce$table$logFC, -log10(qlf.Ca.Vs.Ce$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Cerebellar')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ca.Vs.Ce$table[,"PValue"], breaks=50, main="Caudate vs Cerebellar p-value frequency histogram")

# Caudate vs Cerebellum
qlf.Ca.Vs.Cm <- glmQLFTest(fit, contrast=my.contrasts[,"Ca.Vs.Cm"])
df_Ca.Vs.Cm <- summary(decideTests(qlf.Ca.Vs.Cm))
grid.newpage()
grid.table(df_Ca.Vs.Cm)
plotMD(qlf.At.Vs.Sn)

# Volcano plot
volcanoData <- cbind(qlf.Ca.Vs.Cm$table$logFC, -log10(qlf.Ca.Vs.Cm$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Cerebellum')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ca.Vs.Cm$table[,"PValue"], breaks=50, main="Caudate vs Cerebellum p-value frequency histogram")

# Caudate vs Cortex
qlf.Ca.Vs.Co <- glmQLFTest(fit, contrast=my.contrasts[,"Ca.Vs.Co"])
df_Ca.Vs.Co <- summary(decideTests(qlf.Ca.Vs.Co))
grid.newpage()
grid.table(df_Ca.Vs.Co)
plotMD(qlf.Ca.Vs.Co)

# Volcano plot
volcanoData <- cbind(qlf.Ca.Vs.Co$table$logFC, -log10(qlf.Ca.Vs.Co$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ca.Vs.Co$table[,"PValue"], breaks=50, main="Caudate vs Cortex p-value frequency histogram")

# Caudate vs Frontal Cortex
qlf.Ca.Vs.Fc <- glmQLFTest(fit, contrast=my.contrasts[,"Ca.Vs.Fc"])
df_Ca.Vs.Fc <- summary(decideTests(qlf.Ca.Vs.Fc))
grid.newpage()
grid.table(df_Ca.Vs.Fc)
plotMD(qlf.Ca.Vs.Fc)

# Volcano plot
volcanoData <- cbind(qlf.Ca.Vs.Fc$table$logFC, -log10(qlf.Ca.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Frontal_Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ca.Vs.Fc$table[,"PValue"], breaks=50, main="Caudate vs Frontal Cortex p-value frequency histogram")

# Caudate vs Hipocampus
qlf.Ca.Vs.Hp <- glmQLFTest(fit, contrast=my.contrasts[,"Ca.Vs.Hp"])
df_Ca.Vs.Hp <- summary(decideTests(qlf.Ca.Vs.Hp))
grid.newpage()
grid.table(df_Ca.Vs.Hp)
plotMD(qlf.Ca.Vs.Hp)

# Volcano plot
volcanoData <- cbind(qlf.Ca.Vs.Hp$table$logFC, -log10(qlf.Ca.Vs.Hp[$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ca.Vs.Hp$table[,"PValue"], breaks=50, main="Caudate vs Hippocampus p-value frequency histogram")

# Caudate vs Hypothalamus
qlf.Ca.Vs.Hy <- glmQLFTest(fit, contrast=my.contrasts[,"Ca.Vs.Hy"])
df_Ca.Vs.Hy <- summary(decideTests(qlf.Ca.Vs.Hy))
grid.newpage()
grid.table(df_Ca.Vs.Hy)
plotMD(qlf.Ca.Vs.Hy)

# Volcano plot
volcanoData <- cbind(qlf.Ca.Vs.Hy$table$logFC, -log10(qlf.Ca.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ca.Vs.Hy$table[,"PValue"], breaks=50, main="Caudate vs Hypothalamus p-value frequency histogram")

# Caudate vs Nucleus Accumbens
qlf.Ca.Vs.Na <- glmQLFTest(fit, contrast=my.contrasts[,"Ca.Vs.Na"])
df_Ca.Vs.Na <- summary(decideTests(qlf.Ca.Vs.Na))
grid.newpage()
grid.table(df_Ca.Vs.Na)
plotMD(qlf.Ca.Vs.Na)

# Volcano plot
volcanoData <- cbind(qlf.Ca.Vs.Na$table$logFC, -log10(qlf.Ca.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ca.Vs.Na$table[,"PValue"], breaks=50, main="Caudate vs Nucleus Accumbens p-value frequency histogram")

# Caudate vs Putamen
qlf.Ca.Vs.Pu <- glmQLFTest(fit, contrast=my.contrasts[,"Ca.Vs.Pu"])
df_Ca.Vs.Pu <- summary(decideTests(qlf.Ca.Vs.Pu))
grid.newpage()
grid.table(df_Ca.Vs.Pu)
plotMD(qlf.Ca.Vs.Pu)

# Volcano plot
volcanoData <- cbind(qlf.Ca.Vs.Pu$table$logFC, -log10(qlf.Ca.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ca.Vs.Pu$table[,"PValue"], breaks=50, main="Caudate vs Putamen p-value frequency histogram")

# Caudate vs Spinal Cord
qlf.Ca.Vs.Sp <- glmQLFTest(fit, contrast=my.contrasts[,"Ca.Vs.Sp"])
df_Ca.Vs.Sp <- summary(decideTests(qlf.Ca.Vs.Sp))
grid.newpage()
grid.table(df_Ca.Vs.Sp)
plotMD(qlf.Ca.Vs.Sp)

# Volcano plot
volcanoData <- cbind(qlf.Ca.Vs.Sp$table$logFC, -log10(qlf.Ca.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ca.Vs.Sp$table[,"PValue"], breaks=50, main="Caudate vs Spinal Cord p-value frequency histogram")

# Caudate vs Substantia Nigra
qlf.Ca.Vs.Sn <- glmQLFTest(fit, contrast=my.contrasts[,"Ca.Vs.Sn"])
df_Ca.Vs.Sn <- summary(decideTests(qlf.Ca.Vs.Sn))
grid.newpage()
grid.table(df_Ca.Vs.Sn)
plotMD(qlf.Ca.Vs.Sn)

# Volcano plot
volcanoData <- cbind(qlf.Ca.Vs.Sn$table$logFC, -log10(qlf.Ca.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Caudate-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ca.Vs.Sn$table[,"PValue"], breaks=50, main="Caudate vs Substantia Nigra p-value frequency histogram")

# Cerebellar vs Cerebellum
qlf.Ce.Vs.Cm <- glmQLFTest(fit, contrast=my.contrasts[,"Ce.Vs.Cm"])
df_Ce.Vs.Cm <- summary(decideTests(qlf.Ce.Vs.Cm))
grid.newpage()
grid.table(df_Ce.Vs.Cm)
plotMD(qlf.Ce.Vs.Cm)

# Volcano plot
volcanoData <- cbind(qlf.Ce.Vs.Cm$table$logFC, -log10(qlf.Ce.Vs.Cm$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Cerebellum')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ce.Vs.Cm$table[,"PValue"], breaks=50, main="Cerebellar vs Cerebellum p-value frequency histogram")

# Cerebellar vs Cortex
qlf.Ce.Vs.Co <- glmQLFTest(fit, contrast=my.contrasts[,"Ce.Vs.Co"])
df_Ce.Vs.Co <- summary(decideTests(qlf.Ce.Vs.Co))
grid.newpage()
grid.table(df_Ce.Vs.Co)
plotMD(qlf.Ce.Vs.Co)

# Volcano plot
volcanoData <- cbind(qlf.Ce.Vs.Co$table$logFC, -log10(qlf.Ce.Vs.Co$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ce.Vs.Co$table[,"PValue"], breaks=50, main="Cerebellar vs Cortex p-value frequency histogram")

# Cerebellar vs Frontal Cortex
qlf.Ce.Vs.Fc <- glmQLFTest(fit, contrast=my.contrasts[,"Ce.Vs.Fc"])
df_Ce.Vs.Fc <- summary(decideTests(qlf.Ce.Vs.Fc))
grid.newpage()
grid.table(df_Ce.Vs.Fc)
plotMD(qlf.Ce.Vs.Fc)

# Volcano plot
volcanoData <- cbind(qlf.Ce.Vs.Fc$table$logFC, -log10(qlf.Ce.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Frontal_Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ce.Vs.Fc$table[,"PValue"], breaks=50, main="Cerebellar vs Frontal Cortex p-value frequency histogram")

# Cerebellar vs Hippocampus
qlf.Ce.Vs.Hp <- glmQLFTest(fit, contrast=my.contrasts[,"Ce.Vs.Hp"])
df_Ce.Vs.Hp <- summary(decideTests(qlf.Ce.Vs.Hp))
grid.newpage()
grid.table(df_Ce.Vs.Hp)
plotMD(qlf.Ce.Vs.Hp)

# Volcano plot
volcanoData <- cbind(qlf.Ce.Vs.Hp$table$logFC, -log10(qlf.Ce.Vs.Hp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ce.Vs.Hp$table[,"PValue"], breaks=50, main="Cerebellar vs Hippocampus p-value frequency histogram")

# Cerebellar vs Hypothalamus
qlf.Ce.Vs.Hy <- glmQLFTest(fit, contrast=my.contrasts[,"Ce.Vs.Hy"])
df_Ce.Vs.Hy <- summary(decideTests(qlf.Ce.Vs.Hy))
grid.newpage()
grid.table(df_Ce.Vs.Hy)
plotMD(qlf.Ce.Vs.Hy)

# Volcano plot
volcanoData <- cbind(qlf.Ce.Vs.Hy$table$logFC, -log10(qlf.Ce.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ce.Vs.Hy$table[,"PValue"], breaks=50, main="Cerebellar vs Hypothalamus p-value frequency histogram")

# Cerebellar vs Nucleus Accumbens
qlf.Ce.Vs.Na <- glmQLFTest(fit, contrast=my.contrasts[,"Ce.Vs.Na"])
df_Ce.Vs.Na <- summary(decideTests(qlf.Ce.Vs.Na))
grid.newpage()
grid.table(df_Ce.Vs.Na)
plotMD(qlf.Ce.Vs.Na)

# Volcano plot
volcanoData <- cbind(qlf.Ce.Vs.Na$table$logFC, -log10(qlf.Ce.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ce.Vs.Na$table[,"PValue"], breaks=50, main="Cerebellar vs Nucleus Accumbens p-value frequency histogram")

# Cerebellar vs Putamen
qlf.Ce.Vs.Pu <- glmQLFTest(fit, contrast=my.contrasts[,"Ce.Vs.Pu"])
df_Ce.Vs.Pu <- summary(decideTests(qlf.Ce.Vs.Pu))
grid.newpage()
grid.table(df_Ce.Vs.Pu)
plotMD(qlf.Ce.Vs.Pu)

# Volcano plot
volcanoData <- cbind(qlf.Ce.Vs.Pu$table$logFC, -log10(qlf.Ce.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ce.Vs.Pu$table[,"PValue"], breaks=50, main="Cerebellar vs Putamen p-value frequency histogram")

# Cerebellar vs Spinal Cord
qlf.Ce.Vs.Sp <- glmQLFTest(fit, contrast=my.contrasts[,"Ce.Vs.Sp"])
df_Ce.Vs.Sp <- summary(decideTests(qlf.Ce.Vs.Sp))
grid.newpage()
grid.table(df_Ce.Vs.Sp)
plotMD(qlf.Ce.Vs.Sp)

# Volcano plot
volcanoData <- cbind(qlf.Ce.Vs.Sp$table$logFC, -log10(qlf.Ce.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ce.Vs.Sp$table[,"PValue"], breaks=50, main="Cerebellar vs Spinal Cord p-value frequency histogram")

# Cerebellar vs Substantia Nigra
qlf.Ce.Vs.Sn <- glmQLFTest(fit, contrast=my.contrasts[,"Ce.Vs.Sn"])
df_Ce.Vs.Sn <- summary(decideTests(qlf.Ce.Vs.Sn))
grid.newpage()
grid.table(df_Ce.Vs.Sn)
plotMD(qlf.Ce.Vs.Sn)

# Volcano plot
volcanoData <- cbind(qlf.Ce.Vs.Sn$table$logFC, -log10(qlf.Ce.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cerebellar-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Ce.Vs.Sn$table[,"PValue"], breaks=50, main="Cerebellar vs Substantia Nigra p-value frequency histogram")

# Cortex vs Frontal Cortex
qlf.Co.Vs.Fc <- glmQLFTest(fit, contrast=my.contrasts[,"Co.Vs.Fc"])
df_Co.Vs.Fc <- summary(decideTests(qlf.Co.Vs.Fc))
grid.newpage()
grid.table(df_Co.Vs.Fc)
plotMD(qlf.Co.Vs.Fc)

# Volcano plot
volcanoData <- cbind(qlf.Co.Vs.Fc$table$logFC, -log10(qlf.Co.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortexr-1*Frontal_Cortex')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Co.Vs.Fc$table[,"PValue"], breaks=50, main="Cortex vs Frontal Cortex p-value frequency histogram")

# Cortex vs Hippocampus
qlf.Co.Vs.Hp <- glmQLFTest(fit, contrast=my.contrasts[,"Co.Vs.Hp"])
df_Co.Vs.Hp <- summary(decideTests(qlf.Co.Vs.Hp))
grid.newpage()
grid.table(df_Co.Vs.Hp)
plotMD(qlf.Co.Vs.Hp)

# Volcano plot
volcanoData <- cbind(qlf.Co.Vs.Hp$table$logFC, -log10(qlf.Co.Vs.Hp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Co.Vs.Hp$table[,"PValue"], breaks=50, main="Cortex vs Hippocampus p-value frequency histogram")

# Cortex vs Hypothalamus
qlf.Co.Vs.Hy <- glmQLFTest(fit, contrast=my.contrasts[,"Co.Vs.Hy"])
df_Co.Vs.Hy <- summary(decideTests(qlf.Co.Vs.Hy))
grid.newpage()
grid.table(df_Co.Vs.Hy)
plotMD(qlf.Co.Vs.Hy)

# Volcano plot
volcanoData <- cbind(qlf.Co.Vs.Fc$table$logFC, -log10(qlf.Co.Vs.Fc$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Co.Vs.Hy$table[,"PValue"], breaks=50, main="Cortex vs Hypothalamus p-value frequency histogram")

# Cortex vs Nucleus Accumbens
qlf.Co.Vs.Na <- glmQLFTest(fit, contrast=my.contrasts[,"Co.Vs.Na"])
df_Co.Vs.Na <- summary(decideTests(qlf.Co.Vs.Na))
grid.newpage()
grid.table(df_Co.Vs.Na)
plotMD(qlf.Co.Vs.Na)

# Volcano plot
volcanoData <- cbind(qlf.Co.Vs.Na$table$logFC, -log10(qlf.Co.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Co.Vs.Na$table[,"PValue"], breaks=50, main="Cortex vs Nucleus Accumbens p-value frequency histogram")

# Cortex vs Putamen
qlf.Co.Vs.Pu <- glmQLFTest(fit, contrast=my.contrasts[,"Co.Vs.Pu"])
df_Co.Vs.Pu <- summary(decideTests(qlf.Co.Vs.Pu))
grid.newpage()
grid.table(df_Co.Vs.Pu)
plotMD(qlf.Co.Vs.Pu)

# Volcano plot
volcanoData <- cbind(qlf.Co.Vs.Pu$table$logFC, -log10(qlf.Co.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Co.Vs.Pu$table[,"PValue"], breaks=50, main="Cortex vs Putamen p-value frequency histogram")

# Cortex vs Spinal Cord
qlf.Co.Vs.Sp <- glmQLFTest(fit, contrast=my.contrasts[,"Co.Vs.Sp"])
df_Co.Vs.Sp <- summary(decideTests(qlf.Co.Vs.Sp))
grid.newpage()
grid.table(df_Co.Vs.Sp)
plotMD(qlf.Co.Vs.Sp)

# Volcano plot
volcanoData <- cbind(qlf.Co.Vs.Sp$table$logFC, -log10(qlf.Co.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Co.Vs.Sp$table[,"PValue"], breaks=50, main="Cortex vs Spinal Cord p-value frequency histogram")

# Cortex vs Substantia Nigra
qlf.Co.Vs.Sn <- glmQLFTest(fit, contrast=my.contrasts[,"Co.Vs.Sn"])
df_Co.Vs.Sn <- summary(decideTests(qlf.Co.Vs.Sn))
grid.newpage()
grid.table(df_Co.Vs.Sn)
plotMD(qlf.Co.Vs.Sn)

# Volcano plot
volcanoData <- cbind(qlf.Co.Vs.Sn$table$logFC, -log10(qlf.Co.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Cortex-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Co.Vs.Sn$table[,"PValue"], breaks=50, main="Cortex vs Substantia Nigra p-value frequency histogram")

# Frontal Cortex vs Hippocampus
qlf.Fc.Vs.Hp <- glmQLFTest(fit, contrast=my.contrasts[,"Fc.Vs.Hp"])
df_Fc.Vs.Hp <- summary(decideTests(qlf.Fc.Vs.Hp))
grid.newpage()
grid.table(df_Fc.Vs.Hp)
plotMD(qlf.Fc.Vs.Hp)

# Volcano plot
volcanoData <- cbind(qlf.Fc.Vs.Hp$table$logFC, -log10(qlf.Fc.Vs.Hp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Hippocampus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Fc.Vs.Hp$table[,"PValue"], breaks=50, main="Frontal Cortex vs Hippocampus p-value frequency histogram")

# Frontal Cortex vs Hypothalamus 
qlf.Fc.Vs.Hy <- glmQLFTest(fit, contrast=my.contrasts[,"Fc.Vs.Hy"])
df_Fc.Vs.Hy <- summary(decideTests(qlf.Fc.Vs.Hy))
grid.newpage()
grid.table(df_Fc.Vs.Hy)
plotMD(qlf.Fc.Vs.Hy)

# Volcano plot
volcanoData <- cbind(qlf.Fc.Vs.Hy$table$logFC, -log10(qlf.Fc.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Fc.Vs.Hy$table[,"PValue"], breaks=50, main="Frontal Cortex vs Hypothalamus p-value frequency histogram")

# Frontal Cortex vs Nucleus Accumbens
qlf.Fc.Vs.Na <- glmQLFTest(fit, contrast=my.contrasts[,"Fc.Vs.Na"])
df_Fc.Vs.Na <- summary(decideTests(qlf.Fc.Vs.Na))
grid.newpage()
grid.table(df_Fc.Vs.Na)
plotMD(qlf.Fc.Vs.Na)

# Volcano plot
volcanoData <- cbind(qlf.Fc.Vs.Na$table$logFC, -log10(qlf.Fc.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Fc.Vs.Na$table[,"PValue"], breaks=50, main="Frontal Cortex vs Nucleus Accumbens p-value frequency histogram")

# Frontal Cortex vs Putamen
qlf.Fc.Vs.Pu <- glmQLFTest(fit, contrast=my.contrasts[,"Fc.Vs.Pu"])
df_Fc.Vs.Pu <- summary(decideTests(qlf.Fc.Vs.Pu))
grid.newpage()
grid.table(df_Fc.Vs.Pu)
plotMD(qlf.Fc.Vs.Pu)

# Volcano plot
volcanoData <- cbind(qlf.Fc.Vs.Pu$table$logFC, -log10(qlf.Fc.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Fc.Vs.Pu$table[,"PValue"], breaks=50, main="Frontal Cortex vs Putamen p-value frequency histogram")

# Frontal Corte vs Spinal Cord
qlf.Fc.Vs.Sp <- glmQLFTest(fit, contrast=my.contrasts[,"Fc.Vs.Sp"])
df_Fc.Vs.Sp <- summary(decideTests(qlf.Fc.Vs.Sp))
grid.newpage()
grid.table(df_Fc.Vs.Sp)
plotMD(qlf.Fc.Vs.Sp)

# Volcano plot
volcanoData <- cbind(qlf.Fc.Vs.Sp$table$logFC, -log10(qlf.Fc.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Fc.Vs.Sp$table[,"PValue"], breaks=50, main="Frontal Cortex vs Spinal Cord p-value frequency histogram")

# Frontal Cortex vs Substantia Nigra
qlf.Fc.Vs.Sn <- glmQLFTest(fit, contrast=my.contrasts[,"Fc.Vs.Sn"])
df_Fc.Vs.Sn <- summary(decideTests(qlf.Fc.Vs.Sn))
grid.newpage()
grid.table(df_Fc.Vs.Sn)
plotMD(qlf.Fc.Vs.Sn)

# Volcano plot
volcanoData <- cbind(qlf.Fc.Vs.Sn$table$logFC, -log10(qlf.Fc.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Frontal_Cortex-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Fc.Vs.Sn$table[,"PValue"], breaks=50, main="Frontal Cortex vs Substantia Nigra p-value frequency histogram")

# Hippocampus vs Hypothalamus
qlf.Hp.Vs.Hy <- glmQLFTest(fit, contrast=my.contrasts[,"Hp.Vs.Hy"])
df_Hp.Vs.Hy <- summary(decideTests(qlf.Hp.Vs.Hy))
grid.newpage()
grid.table(df_Hp.Vs.Hy)
plotMD(qlf.Hp.Vs.Hy)

# Volcano plot
volcanoData <- cbind(qlf.Hp.Vs.Hy$table$logFC, -log10(qlf.Hp.Vs.Hy$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus-1*Hypothalamus')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Hp.Vs.Hy$table[,"PValue"], breaks=50, main="Hippocampus vs Hypothalamus p-value frequency histogram")

# Hippocampus vs Nucleus Accumbens
qlf.Hp.Vs.Na <- glmQLFTest(fit, contrast=my.contrasts[,"Hp.Vs.Na"])
df_Hp.Vs.Na <- summary(decideTests(qlf.Hp.Vs.Na))
grid.newpage()
grid.table(df_Hp.Vs.Na)
plotMD(qlf.Hp.Vs.Na)

# Volcano plot
volcanoData <- cbind(qlf.Hp.Vs.Na$table$logFC, -log10(qlf.Hp.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Hp.Vs.Na$table[,"PValue"], breaks=50, main="Hippocampus vs Nucleus Accumbens p-value frequency histogram")

# Hippocampus vs Putamen
qlf.Hp.Vs.Pu <- glmQLFTest(fit, contrast=my.contrasts[,"Hp.Vs.Pu"])
df_Hp.Vs.Pu <- summary(decideTests(qlf.Hp.Vs.Pu))
grid.newpage()
grid.table(df_Hp.Vs.Pu)
plotMD(qlf.Hp.Vs.Pu)

# Volcano plot
volcanoData <- cbind(qlf.Hp.Vs.Pu$table$logFC, -log10(qlf.Hp.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Hp.Vs.Pu$table[,"PValue"], breaks=50, main="Hippocampus vs Putamen p-value frequency histogram")

# Hippocampus vs Spinal Cord
qlf.Hp.Vs.Sp <- glmQLFTest(fit, contrast=my.contrasts[,"Hp.Vs.Sp"])
df_Hp.Vs.Sp <- summary(decideTests(qlf.Hp.Vs.Sp))
grid.newpage()
grid.table(df_Hp.Vs.Sp)
plotMD(qlf.Hp.Vs.Sp)

# Volcano plot
volcanoData <- cbind(qlf.Hp.Vs.Sp$table$logFC, -log10(qlf.Hp.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Hp.Vs.Sp$table[,"PValue"], breaks=50, main="Hippocampus vs Spinal Cord p-value frequency histogram")

# Hippocampus vs Substantia Nigra 
qlf.Fc.Vs.Sn <- glmQLFTest(fit, contrast=my.contrasts[,"Fc.Vs.Sn"])
df_Hp.Vs.Sn <- summary(decideTests(qlf.Hp.Vs.Sn))
grid.newpage()
grid.table(df_Hp.Vs.Sn)
plotMD(qlf.Hp.Vs.Sn)

# Volcano plot
volcanoData <- cbind(qlf.Hp.Vs.Sn$table$logFC, -log10(qlf.Hp.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hippocampus-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Hp.Vs.Sn$table[,"PValue"], breaks=50, main="Hippocampus vs Substantia Nigra p-value frequency histogram")

# Hypothalamus vs Nucleus Accumbens
qlf.Hy.Vs.Na <- glmQLFTest(fit, contrast=my.contrasts[,"Hy.Vs.Na"])
df_Hy.Vs.Na <- summary(decideTests(qlf.Hy.Vs.Na))
grid.newpage()
grid.table(df_Hy.Vs.Na)
plotMD(qlf.Hy.Vs.Na)

# Volcano plot
volcanoData <- cbind(qlf.Hy.Vs.Na$table$logFC, -log10(qlf.Hy.Vs.Na$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hypothalamus-1*Nucleus_Accumbens')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Hy.Vs.Na$table[,"PValue"], breaks=50, main="Hypothalamus vs Nucleus Accumbens p-value frequency histogram")

# Hypothalamus vs Putamen
qlf.Hy.Vs.Pu <- glmQLFTest(fit, contrast=my.contrasts[,"Hy.Vs.Pu"])
df_Hy.Vs.Pu <- summary(decideTests(qlf.Hy.Vs.Pu))
grid.newpage()
grid.table(df_Hy.Vs.Pu)
plotMD(qlf.Hy.Vs.Pu)

# Volcano plot
volcanoData <- cbind(qlf.Hy.Vs.Pu$table$logFC, -log10(qlf.Hy.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hypothalamus-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Hy.Vs.Pu$table[,"PValue"], breaks=50, main="Hypothalamus vs Putamen p-value frequency histogram")

# Hypothalamus vs Spinal Cord
qlf.Hy.Vs.Sp <- glmQLFTest(fit, contrast=my.contrasts[,"Hy.Vs.Sp"])
df_Hy.Vs.Sp <- summary(decideTests(qlf.Hy.Vs.Sp))
grid.newpage()
grid.table(df_Hy.Vs.Sp)
plotMD(qlf.Hy.Vs.Sp)

# Volcano plot
volcanoData <- cbind(qlf.Hy.Vs.Sp$table$logFC, -log10(qlf.Hy.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hypothalamus-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Hy.Vs.Sp$table[,"PValue"], breaks=50, main="Hypothalamus vs Spinal Cord p-value frequency histogram")

# Hypothalamus Substantia Nigra
qlf.Hy.Vs.Sn <- glmQLFTest(fit, contrast=my.contrasts[,"Hy.Vs.Sn"])
df_Hy.Vs.Sn <- summary(decideTests(qlf.Hy.Vs.Sn))
grid.newpage()
grid.table(df_Hy.Vs.Sn)
plotMD(qlf.Hy.Vs.Sn)

# Volcano plot
volcanoData <- cbind(qlf.Hy.Vs.Sn$table$logFC, -log10(qlf.Hy.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Hypothalamus-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Hy.Vs.Sn$table[,"PValue"], breaks=50, main="Hypothalamus vs Substantia Nigra p-value frequency histogram")

# Nucleus Accumbens vs Putamen
qlf.Na.Vs.Pu <- glmQLFTest(fit, contrast=my.contrasts[,"Na.Vs.Pu"])
df_Na.Vs.Pu <- summary(decideTests(qlf.Na.Vs.Pu))
grid.newpage()
grid.table(df_Na.Vs.Pu)
plotMD(qlf.Na.Vs.Pu)

# Volcano plot
volcanoData <- cbind(qlf.Na.Vs.Pu$table$logFC, -log10(qlf.Na.Vs.Pu$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Nucleus_Accumbens-1*Putamen')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Na.Vs.Pu$table[,"PValue"], breaks=50, main="Nucleus Accumbens vs Putamen p-value frequency histogram")

# Nucleus Accumbens vs Spinal Cord
qlf.Na.Vs.Sp <- glmQLFTest(fit, contrast=my.contrasts[,"Na.Vs.Sp"])
df_Na.Vs.Sp <- summary(decideTests(qlf.Na.Vs.Sp))
grid.newpage()
grid.table(df_Na.Vs.Sp)
plotMD(qlf.Na.Vs.Sp)

# Volcano plot
volcanoData <- cbind(qlf.Na.Vs.Sp$table$logFC, -log10(qlf.Na.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Nucleus_Accumbens-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Na.Vs.Sp$table[,"PValue"], breaks=50, main="Nucleus Accumbens vs Spinal Cord p-value frequency histogram")

# Nucleus Accumbens vs Substantia Nigra
qlf.Na.Vs.Sn <- glmQLFTest(fit, contrast=my.contrasts[,"Na.Vs.Sn"])
df_Na.Vs.Sn <- summary(decideTests(qlf.Na.Vs.Sn))
grid.newpage()
grid.table(df_na.Vs.Sn)
plotMD(qlf.Na.Vs.Sn)

# Volcano plot
volcanoData <- cbind(qlf.Na.Vs.Sn$table$logFC, -log10(qlf.Na.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Nucleus_Accumbens-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Na.Vs.Sn$table[,"PValue"], breaks=50, main="Nucleus Accumbens vs Substantia Nigra p-value frequency histogram")

# Putamen vs Spinal Cord
qlf.Pu.Vs.Sp <- glmQLFTest(fit, contrast=my.contrasts[,"Pu.Vs.Sp"])
df_Pu.Vs.Sp <- summary(decideTests(qlf.Pu.Vs.Sp))
grid.newpage()
grid.table(df_Pu.Vs.Sp)
plotMD(qlf.Pu.Vs.Sp)

# Volcano plot
volcanoData <- cbind(qlf.Pu.Vs.Sp$table$logFC, -log10(qlf.Pu.Vs.Sp$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Putamen-1*Spinal_Cord')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Pu.Vs.Sp$table[,"PValue"], breaks=50, main="Putamen vs Spinal Cord p-value frequency histogram")

# Putamen vs Substantia Nigra
qlf.Pu.Vs.Sn <- glmQLFTest(fit, contrast=my.contrasts[,"Pu.Vs.Sn"])
df_Pu.Vs.Sn <- summary(decideTests(qlf.Pu.Vs.Sn))
grid.newpage()
grid.table(df_Pu.Vs.Sn)
plotMD(qlf.Pu.Vs.Sn)

# Volcano plot
volcanoData <- cbind(qlf.Pu.Vs.Sn$table$logFC, -log10(qlf.Pu.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Putamen-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Pu.Vs.Sn$table[,"PValue"], breaks=50, main="Putamen vs Substantia Nigra p-value frequency histogram")

# Spinal Cord vs Substantia Nigra
qlf.Sp.Vs.Sn <- glmQLFTest(fit, contrast=my.contrasts[,"Sp.Vs.Sn"])
df_Sp.Vs.Sn <- summary(decideTests(qlf.Sp.Vs.Sn))
grid.newpage()
grid.table(df_Sp.Vs.Sn)
plotMD(qlf.Sp.Vs.Sn)

# Volcano plot
volcanoData <- cbind(qlf.Sp.Vs.Sn$table$logFC, -log10(qlf.Sp.Vs.Sn$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Spinal_Cord-1*Substantia_Nigra')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Sp.Vs.Sn$table[,"PValue"], breaks=50, main="Spinal Cord vs Substantia Nigra p-value frequency histogram")

dev.off()

