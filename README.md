# GTEx_Brains

Quantify differential gene expression in RNA Seq data from the Genotype Tissue Expression (GTEx) project between individuals, sexes, and tissues in 11 brain regions using various DGX pipelines.

Regions: Anterior, Amygdala, Caudate, Cerebellum, Cortex, Hippocampus, Hypothalamus, Nucleus Accumbens, Putamen, Spinal Cord, Substantia Nigra.

Sample Replicates: Cerebellum/Cortex samples are replicates of Cerebellar/Frontal Cortex samples, respectively. Cortex and Cerebellum samples differ from other samples in time of collection and manner in which they were preserved. 

Quality Control: FastQC, MultiQC, Trimmomatic, Trimmed-FastQC, Trimmed-MultiQC

Quantification: Salmon, Hisat2

Differential Expression: DESeq2, EdgeR, Limma/Voom, IsoformSwitchAnalyzeR
