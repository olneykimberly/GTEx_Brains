# Write table to see what types of tissue samples are available per individuals, including whole blood samples.

library(jsonlite)
library(stringr)

# Read in Metadata.csv
meta <- read.csv(file.path("/scratch/mjpete11/GTEx/", "No_Dup_Metadata.csv"), header = TRUE)

# Read in config containing dictionary mapping GTEx IDs to Illumina SRA IDs 
map <- read_json(path="/scratch/mjpete11/GTEx/Spinal_Cord/Male_Spinal_Cord.config.json")

# Convert config to list of lists
map_lst <- fromJSON("/scratch/mjpete11/GTEx/Spinal_Cord/Male_Spinal_Cord.config.json")

# Read in list of whole blood SRA IDs
SRA_Whole_Blood <- list.files(path = "/data/storage/public/dbgap-8834/whole_blood_rna/whole_blood", pattern = "*_1.fastq", all.files = TRUE)

# Keep only run ID
SRA_Whole_Blood <- str_remove_all(SRA_Whole_Blood, "_1.fastq")

# Make named character vector: Names are GTEx IDs, values are SRA IDs
All_GTEx_IDs <- sapply(map_lst, '[[', 9)

# Check that the whole blood samples on agave are in the GTEx sample list
SRA_Whole_Blood %in% All_GTEx_IDs # TRUE

# Get whole blood GTEx IDs
GTEx_Whole_Blood <- names(All_GTEx_IDs[(which(SRA_Whole_Blood %in% All_GTEx_IDs))])
 
# Write table
tmp_mat <- rbind(meta[,1:2],cbind(Sample=GTEx_Whole_Blood,Tissue="Blood"))
tmp_mat$Sample <- str_extract(tmp_mat$Sample,"^GTEX-[0-9A-Z]{4}")
tmp_mat <- unique(tmp_mat)
sample_table <- as.data.frame.matrix(table(tmp_mat))
write.csv(sample_table,file = "Blood_Brain_Power.csv")

