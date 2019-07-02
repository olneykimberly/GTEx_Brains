# Wilcox distribution of age of the samples

library(dplyr)
library(ggplot2)

# Constants
PLOT_DIR <- "/scratch/mjpete11/GTEx/Data_Exploration/Hist/"
PLOT_NAME <- "Age_Histograms.pdf"

# Read Metadata CSV.                                                            
Samples = read.csv(file.path("/scratch/mjpete11/GTEx/Metadata/", "Metadata.csv"), header = TRUE)
Samples                                                                         

# Vector of ages of female and male samples
Females <- Samples %>%
             select(Age, Sex) %>%
               filter(Sex == "Female") %>%
                 unlist()

Males <- Samples %>%
           select(Age, Sex) %>%
             filter(Sex == "Male") %>%
               unlist()

# Wilcox test
# H0: There is no difference in the age distribution of female and male samples
w_test <- wilcox.test(Females, Males, conf.level= 0.95)

w_test$p.value

# Result for config without sample duplicates: 2.403928e-23

# Overlaid density plots
ggplot(Samples, aes(x=Samples$Age, fill=Samples$Sex)) +
  geom_density(alpha=.3) + ggtitle("Density Plot of Sample Age Distribution") +
  xlab("Age") + ylab("Density") + scale_fill_manual(name = "Sex", values=c("blue", "green"))

# Organized samples by tissue type into list of dfs
Meta <- list()
for(i in 1:length(levels(Samples$Tissue))){
  Meta[[i]] <- Samples[Samples$Tissue == levels(Samples$Tissue)[i],]
}
names(Meta) <- levels(Samples$Tissue)

# X lims
max(Samples$Age) # 70
min(Samples$Age) # 20

# Histogram plots
setwd(PLOT_DIR)
pdf(PLOT_NAME)

x_axis_labels <- seq(min(Samples[,'Age']), max(Samples[,'Age']), 5)
y_axis_labels <- seq(0, 30, 1)

Hist_Func <- function(a, b){
  Hist_Plots <- ggplot(a, aes(x=Age, fill=Sex)) +
    geom_histogram(data=subset(a, Sex=='Female'), aes(fill=Sex), alpha=0.4, binwidth=1, colour='gray50') +
    geom_histogram(data=subset(a, Sex=='Male'), aes(fill=Sex), alpha=0.4, binwidth=1, colour='gray50') + 
    scale_x_continuous(labels=x_axis_labels, breaks=x_axis_labels) +
    scale_y_continuous(labels=y_axis_labels, breaks=y_axis_labels) +
    scale_fill_manual(name="Sex", values=c("blue", "green"), labels=c("Female", "Male")) +
    ggtitle(paste(b, "Sample Age Histogram")) +
    xlab("Age") + ylab("Count") 
}

Map(Hist_Func, a = Meta, b = names(Meta))

dev.off()


