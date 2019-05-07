# This script is for plotting PCA dimensions 1:4 for the top 1000 most variable gene per tissue.

METADATA = file.path("/scratch/mjpete11/GTEx/", "Metadata.csv")
COUNTS = file.path("/scratch/mjpete11/GTEx/Data_Exploration/", "Count_Matrix.tsv")
PLOT_DIR = "PCA_Tissue_Plots/"

# Load packages                                                                 
library(readr)
library(stringr)
library(ggplot2)
library(dplyr)

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

# Filter for most variable genes
var_genes <- apply(cts, 1, var) # 1 is to  set margin to rows

# Select 1000 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE)[1:1000]) 

# Subset the most variable genes
highly_variable_cpm <- cts[select_var,]
dim(highly_variable_cpm)

# Create list of data frames containing counts for each tissue
var_count <- list()
for(i in 1:length(tissue_lst)){
  var_count[[i]] <- highly_variable_cpm[,which(colnames(highly_variable_cpm) %in% tissue_lst[[i]])]
}

# Swap rows and columns in each df in list
var_count <- lapply(var_count, t)

# Rename lists in list
names(var_count) <- paste("count_",names(tissue_lst),sep="")

# Euclidean distances between the rows
count_dist <- lapply(var_count, dist)

# Calculate PCA
fitD <- lapply(count_dist, function(y){
  cmdscale(y, eig=T, k=4) # Set number of dimensions
})

# Add metadata to PCA data
# Get metadata for each tissue and store dataframes in list
meta <- list()
for(i in 1:length(levels(samples$Tissue))){
  meta[[i]] <- samples[samples$Tissue == levels(samples$Tissue)[i],]
}
names(meta) <- levels(samples$Tissue)

# Merge metadata and PCA data into single dataframe
# e.g. combine two lists of dataframes into a list of merged pairs of dfs
mds <- list()
for(i in 1:length(fitD)){
  # make data frame with samples as a column
  tmpA <- data.frame(Sample=rownames(fitD[[i]][[1]]),fitD[[i]][[1]])
  # save meta as a tmp dataframe
  tmpB <- meta[[i]]
  # refactor tmpB data frame to avoid errors
  tmpB$Sample <- factor(tmpB$Sample)
  # join tables together by samples name
  mds[[i]] <- left_join(tmpA,tmpB)
  # match unnamed columns
  tmpC <- str_which(colnames(mds[[i]]),pattern = "^[A-Z][0-9]+")
  # rename columns
  colnames(mds[[i]])[tmpC] <- paste("Dim_",1:length(tmpC),sep = "")
}
names(mds) <- names(meta)


# PCA plots
setwd(PLOT_DIR)

# Make plots and store in list
plot_d2_sex_tissue <- list()
plot_d2_age_tissue <- list()
plot_d4_sex_tissue <- list()
plot_d4_age_tissue <- list()

 for (i in 1:length(mds)) {
   plot_d2_sex_tissue[[i]] <- ggplot(mds[[i]], aes(x=Dim_1, y=Dim_2, shape=Sex, color=Sex, label=Tissue)) +
                              geom_point(size=2, aes(fill = Sex)) + ggtitle(paste("PCA plot",names(mds)[i],sep = " : " )) +
                              scale_shape_manual(values = c(21,22)) + scale_color_manual(values=c("blue3", "green3")) +
                              scale_fill_manual(values=c("blue3", "green3"))
   plot_d2_age_tissue[[i]] <- ggplot(mds[[i]], aes(x=Dim_1, y=Dim_2, color=Age)) +
                              geom_point(size=2) + ggtitle(paste("PCA plot",names(mds)[i],sep = " : " )) + 
                              scale_color_gradient(low="black", high="lightblue")
   plot_d4_sex_tissue[[i]] <- ggplot(mds[[i]], aes(x=Dim_3, y=Dim_4, shape=Sex, color=Sex, label=Tissue)) +
                              geom_point(size=2, aes(fill = Sex)) + ggtitle(paste("PCA plot",names(mds)[i],sep = " : " )) +
                              scale_shape_manual(values = c(21,22)) + scale_color_manual(values=c("blue3", "green3")) +
                              scale_fill_manual(values=c("blue3", "green3"))
   plot_d4_age_tissue[[i]] <- ggplot(mds[[i]], aes(x=Dim_3, y=Dim_4, color=Age)) +
                              geom_point(size=2) + ggtitle(paste("PCA plot",names(mds)[i],sep = " : " )) + 
                              scale_color_gradient(low="black", high="lightblue")
 }

 names(plot_d2_sex_tissue) <- names(mds)
 names(plot_d2_age_tissue) <- names(mds)
 names(plot_d4_sex_tissue) <- names(mds)
 names(plot_d4_age_tissue) <- names(mds)
 
# Print plots of  dimensions 1 and 2: Sex and Tissue

for (i in 1:length(plot_d2_sex_tissue)) {
  file_name = paste(names(mds)[i], "_d2_sex_tissue", ".pdf", sep="")
  pdf(file_name)
  print(plot_d2_sex_tissue[[i]])
  dev.off()
}

# Print plots of dimensions 1 and 2: Age and Tissue
for (i in 1:length(plot_d2_age_tissue)) {
  file_name = paste(names(mds)[i], "_d2_age_tissue", ".pdf", sep="")
  pdf(file_name)
  print(plot_d2_age_tissue[[i]])
  dev.off()
}

# Print plots of dimension 3 and 4: Sex and Tissue
for (i in 1:length(plot_d4_sex_tissue)) {
  file_name = paste(names(mds)[i], "_d4_sex_tissue", ".pdf", sep="")
  pdf(file_name)
  print(plot_d4_sex_tissue[[i]])
  dev.off()
}

# Print plots of dimension 3 and 4: Age and Tissue
for (i in 1:length(plot_d4_age_tissue)) {
  file_name = paste(names(mds)[i], "_d4_age_tissue", ".pdf", sep="")
  pdf(file_name)
  print(plot_d4_age_tissue[[i]])
  dev.off()
}

# Compile all .pdf  files in directory into single file
system("convert *.pdf Combined_Plots_Tissue.pdf")