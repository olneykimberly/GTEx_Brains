# This script looks at differential gene expression between females and males across all tissue.

METADATA = "/scratch/mjpete11/GTEx/Metadata.csv"
COUNTS = "/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrix.tsv"
PLOT_DIR =  "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex/Gene/Normalized/GLM_F_Test/Plots/"

# Load packages                                                                 
library(tximport)                                                               
library(GenomicFeatures)                                                        
library(edgeR) 
library(stringr)
library(gridExtra)
library(grid)

# Read Metadata CSV.                                                            
samples = read.csv(METADATA, header = TRUE)
samples

# Set rownames of metadata object equal to sample names.                        
rownames(samples) <- samples$Sample                                             

# Read in counts 
cts <- read.csv(COUNTS, sep = "\t")

# Replace . to - in colnames
colnames(cts) <- str_replace_all(colnames(cts),pattern = "\\.","-")

# Create design matrix: Step 1
# 26 combos of tissue and sex
Group <- factor(paste(samples$Tissue, samples$Sex, sep="."))
cbind(samples,group=Group)

# Create DGEList object: Sex only
y <- DGEList(cts, group=samples$Sex)

#Create deisgn matrix: Step 2
design <- model.matrix(~0+group, data=y$samples) # No intercept
colnames(design) <- levels(y$samples$group)
design

# Filter out lowly expressed genes. (< 7 counts)
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]

# TMM is the recommended for most RNA-Seq data where the majority (more than half) of the genes 
# are believed not differentially expressed between any pair of the samples.
y <- calcNormFactors(y)
y$samples

# Estimate common dispersion and tagwise dispersions in one run (recommended)
# The square root of the common dispersion gives the coefficient of variation of biological variation.
y <- estimateDisp(y, design, robust=TRUE)
sqrt(y$common.dispersion) 

# Plot
setwd(PLOT_DIR)

pdf('GLM_F_Test_Fe.Vs.Ma.pdf')

# Plots
plotBCV(y, main="Biological Coefficient of Variation Plot")

# Test for DGX genes using quasi-liklihood method.
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit, main="Quasi-Likelihood Disperson")

# Make contrasts: Female vs Male
my.contrasts <- makeContrasts(Fe.Vs.Ma = Female - Male, 
                              Ma.Vs.Fe = Male - Female, levels = design)

# Test for DGX genes using quasi-liklihood method.
# Female vs Male
qlf.Fe.Vs.Ma <- glmQLFTest(fit, contrast=my.contrasts[,"Fe.Vs.Ma"])
df <- summary(decideTests(qlf.Fe.Vs.Ma))
grid.newpage() # To keep summary table from plotting on top of other plots
grid.table(df)
plotMD(qlf.Fe.Vs.Ma)

# Volcano plot
volcanoData <- cbind(qlf.Fe.Vs.Ma$table$logFC, -log10(qlf.Fe.Vs.Ma$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Female-1*Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(qlf.Fe.Vs.Ma$table[,"PValue"], breaks=50, main="Female vs Male p-value frequency histogram")

dev.off()
