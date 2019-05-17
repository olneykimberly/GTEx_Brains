# Create config of samples that worked with Salmon, without sample replicates and only sample replicates.

library(jsonlite)
library(stringr)

Quant_Files <- read_json("/scratch/mjpete11/GTEx/Dup_Quantification.config.json")

# Check if samples have same individual ID
# If same ID, keep only one.
No_Reps <- lapply(Quant_Files, function(y){
    c <- y[-which(duplicated(str_extract(y,"GTEX-[0-9A-Z]+")))]
    if (length(c)==0) {y} else {c}
})

# Write config
write_json(No_Reps,"/scratch/mjpete11/GTEx/Quantification.config.json")

# Identify sample replicates
Filt_Quant <- Quant_Files[-c(1:3)];

Reps <- lapply(Filt_Quant, function(y){
  c <- y[which(duplicated(str_extract(y,"GTEX-[0-9A-Z]+")))]
  if (length(c)==0) {NA} else {c}
})

# Write config
write_json(Reps,"/scratch/mjpete11/GTEx/Rep_Only_Quantification.config.json")


