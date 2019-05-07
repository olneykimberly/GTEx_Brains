# This script looks at differential gene expression between females and males across all tissue.

METADATA = "/scratch/mjpete11/GTEx/Metadata.csv"
COUNTS = "/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrix.tsv"
PLOT_DIR =  "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex/Gene/Unnormalized/GLM_Ratio_Test/Plots/"

# Load packages                                                                 
library(tximport)                                                               
library(GenomicFeatures)                                                        
library(edgeR) 
library(stringr)
library(gridExtra)
library(grid)

# Read Metadata CSV.                                                            
samples = read.csv(METADATA, header = TRUE)

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

# Filter out lowly expressed genes. (< 7 counts)
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]

# Estimate common dispersion and tagwise dispersions in one run (recommended)
# The square root of the common dispersion gives the coefficient of variation of biological variation.
y <- estimateDisp(y, design, robust=TRUE)

# Test for DGX using liklihood ratio test.
fit <- glmFit(y, design, robust=TRUE)

# Make contrasts: Female vs Male
my.contrasts <- makeContrasts(Fe.Vs.Ma = Female - Male,
                              Ma.Vs.Fe = Male - Female, levels = design)

# Plot
setwd(PLOT_DIR)

pdf('GLM_Ratio_Test_Fe.Vs.Ma.pdf')

# Test for DGX genes using ratio test.
# Female vs Male
lrt.Fe.Vs.Ma <- glmLRT(fit, contrast=my.contrasts[,"Fe.Vs.Ma"])
df <- summary(decideTests(lrt.Fe.Vs.Ma))
grid.newpage() 
grid.table(df)
plotMD(lrt.Fe.Vs.Ma)

# Volcano plot
volcanoData <- cbind(lrt.Fe.Vs.Ma$table$logFC, -log10(lrt.Fe.Vs.Ma$table[,"PValue"]))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, main='Volcano plot: 1*Female-1*Male')
abline(a=1.30102999566, b=0, col="blue") # Set intercept equal to p = -log10(0.05)
abline(v=0, col="red")

# Plot unadjusted p-values
hist(lrt.Fe.Vs.Ma$table[,"PValue"], breaks=50, main="Female vs Male p-value frequency histogram")

dev.off()
