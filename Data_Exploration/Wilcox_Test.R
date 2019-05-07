# Wilcox distribution of age of the samples

library(dplyr)
library(ggplot2)

# Read Metadata CSV.                                                            
samples = read.csv(file.path("/scratch/mjpete11/GTEx/", "No_Dup_Metadata.csv"), header = TRUE)
samples                                                                         

# Vector of ages of female and male samples
females <- samples %>%
             select(Age, Sex) %>%
               filter(Sex == "Female") %>%
                 unlist()

males <- samples %>%
           select(Age, Sex) %>%
             filter(Sex == "Male") %>%
               unlist()

# Wilcox test
# H0: There is no difference in the age distribution of female and male samples
w_test <- wilcox.test(females, males, conf.level= 0.95)

w_test$p.value

# Result for config with sample duplicates: 1.909081e-20 
# Result for config without sample duplicates: 2.403928e-23

# Overlaid density plots
ggplot(samples, aes(x=samples$Age, fill=samples$Sex)) +
  geom_density(alpha=.3) + ggtitle("Density Plot of Sample Age Distribution") +
  xlab("Age") + ylab("Density") + scale_fill_manual(name = "Sex", values=c("blue", "green"))



