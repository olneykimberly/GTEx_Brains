# Create config of samples that worked with Salmon, without sample replicates.

library(jsonlite)
library(stringr)

x <- read_json("/scratch/mjpete11/GTEx/Dup_Quantification.config.json")

# Check if samples have same individual ID
# If same ID, keep only one.
tmp <- lapply(x,function(y){
    c <- y[-which(duplicated(str_extract(y,"GTEX-[0-9A-Z]+")))]
    if(length(c)==0){y}else{c}
})

write_json(tmp,"/scratch/mjpete11/GTEx/No_Dup_Quantification.config.json")


