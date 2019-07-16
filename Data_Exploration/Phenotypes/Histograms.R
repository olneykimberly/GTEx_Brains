# Histograms of sample phenotypes labelled by sex

library(dplyr)
library(ggplot2)

# Constants
PLOT_DIR <- "/scratch/mjpete11/GTEx/Data_Exploration/Phenotypes/"
AGE_HIST <- "Age_Histograms.pdf"
RACE_HIST <- "Race_Histograms.pdf"
ETH_HIST <- "Ethnicity_Histograms.pdf"
AGE_MATCHED_HIST <- "Age_Matched_Histograms.pdf"
AGE_DENS <- "Age_Density_Plots.pdf"
AGE_MATCHED_DENS <- "Age_Matched_Density_Plots.pdf"


# Read sample metadata.                                                            
Samples <- read.csv(file.path("/scratch/mjpete11/GTEx/Metadata/", "Metadata.csv"), header = TRUE)
Age_Matched <- read.csv(file.path("/scratch/mjpete11/GTEx/Metadata/", "Age_Matched_Metadata.csv"), header = TRUE)

# Min and max of sample ages
min(Samples$Age) # 20
max(Samples$Age) # 70

min(Age_Matched$Age) # 51
max(Age_Matched$Age) # 70

# Overlaid density plots for all samples
ggplot(Samples, aes(x=Samples$Age, fill=Samples$Sex)) +
 geom_density(alpha=.3) + ggtitle("Density Plot of Sample Age Distribution") +
 xlab("Age") + ylab("Density") + scale_fill_manual(name = "Sex", values=c("blue", "green"))

# Same; age-matched
ggplot(Age_Matched, aes(x=Age_Matched$Age, fill=Age_Matched$Sex)) +
  geom_density(alpha=.3) + ggtitle("Density Plot of Age-Matched Sample Age Distribution") +
  xlab("Age") + ylab("Density") + scale_fill_manual(name = "Sex", values=c("blue", "green"))

# Organize samples by tissue type into list of dfs
Meta_Samples <- list()
for(i in 1:length(levels(Samples$Tissue))){
  Meta_Samples[[i]] <- Samples[Samples$Tissue == levels(Samples$Tissue)[i],]
}
names(Meta_Samples) <- levels(Samples$Tissue)

# Same for age-matched samples
Meta_Age_Matched <- list()
for(i in 1:length(levels(Age_Matched$Tissue))){
  Meta_Age_Matched[[i]] <- Age_Matched[Age_Matched$Tissue == levels(Age_Matched$Tissue)[i],]
}
names(Meta_Age_Matched) <- levels(Age_Matched$Tissue)

# Plots
setwd(PLOT_DIR)

# Function to plot histograms on top of each other
x_axis_labels <- seq(min(Samples[,'Age']), max(Samples[,'Age']), 5)
y_axis_labels <- seq(0, 30, 1)

Hist_Func <- function(META, NAMES, TITLE){
  Hist_Plots <- ggplot(META, aes(x=Age, fill=Sex)) +
    geom_histogram(data=subset(META, Sex=='Female'), aes(fill=Sex), alpha=0.4, binwidth=1, colour='gray50') +
    geom_histogram(data=subset(META, Sex=='Male'), aes(fill=Sex), alpha=0.4, binwidth=1, colour='gray50') + 
    scale_x_continuous(labels=x_axis_labels, breaks=x_axis_labels) +
    scale_y_continuous(labels=y_axis_labels, breaks=y_axis_labels) +
    scale_fill_manual(name="Sex", values=c("blue", "green"), labels=c("Female", "Male")) +
    ggtitle(paste(NAMES, TITLE)) +
    xlab("Age") + ylab("Count") 
}

# Plot sample histograms
pdf(AGE_HIST)
Map(Hist_Func, META = Meta_Samples, NAMES = names(Meta_Samples), TITLE = "Sample Age Histogram")
dev.off()

# Plot age-matched sample histogram
pdf(AGE_MATCHED_HIST)
Map(Hist_Func, META = Meta_Age_Matched, NAMES = names(Meta_Age_Matched), TITLE = "Age-Matched Sample Histogram")
dev.off()

# Density plots by tissue
Density_Func <- function(META, NAMES, TITLE){
  Density_Plots <- ggplot(META, aes(x=Age, fill=Sex)) +
    geom_density(data=subset(META), alpha=.3) + 
    scale_fill_manual(name = "Sex", values=c("blue", "green")) +  
    ggtitle(paste(NAMES, TITLE)) +
    xlab("Age") + ylab("Density") 
}

# Plot samples
pdf(AGE_DENS)
Map(Density_Func, META = Meta_Samples, NAMES = names(Meta_Samples), TITLE = "Sample Age Density Plot")
dev.off()

# Plot age-matched
pdf(AGE_MATCHED_DENS)
Map(Density_Func, META = Meta_Age_Matched, NAMES = names(Meta_Age_Matched), TITLE = "Age-Matched Sample Density Plot")
dev.off()

# Race histogram plots
y_axis_labels <- seq(0, 100, 2)

Hist_Func <- function(a, b){
  a$Race<-as.character(a$Race)
  Hist_Plots <-ggplot(a, aes(x=Race, fill=Sex)) +
    geom_bar(stat = "count", position = "dodge") +
    scale_fill_manual(name="Sex", values=c("blue", "green"), labels=c("Female", "Male")) +
    scale_y_continuous(labels=y_axis_labels, breaks=y_axis_labels) +
    ggtitle(paste(b, "Sample Race Histogram")) +
    xlab("Race") + ylab("Count")
}

pdf(RACE_HIST)
Map(Hist_Func, a = Meta, b = names(Meta)) 
dev.off()

# Ethnicity histogram plots
y_axis_labels <- seq(0, 100, 2) 

Hist_Func <- function(a, b){
  a$Ethnicity<-as.character(a$Ethnicity)
  Hist_Plots <-ggplot(a, aes(x=Ethnicity, fill=Sex)) +
    geom_bar(stat = "count", position = "dodge") +
    scale_fill_manual(name="Sex", values=c("blue", "green"), labels=c("Female", "Male")) +
    scale_y_continuous(labels=y_axis_labels, breaks=y_axis_labels) +
    ggtitle(paste(b, "Sample Ethnicity Histogram")) +
    xlab("Ethnicity") + ylab("Count")
}

pdf(ETH_HIST)
Map(Hist_Func, a = Meta, b = names(Meta))
dev.off()

