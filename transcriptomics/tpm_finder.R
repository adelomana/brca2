### if you need to install libraries, use these lines
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# BiocManager::install("tximport")
# BiocManager::install("biomaRt")
# BiocManager::install("rhdf5")
#####################################################

###
### 0. load libraries
###
library(DESeq2)
library(tximport)
library(biomaRt)

###
### 1. user-defined variables
###
results_dir = '/home/adrian/projects/brca2/results/tpm/'
kallisto_dir = "/home/adrian/projects/brca2/results/kallisto/kallisto.100"
metadata_file = '/home/adrian/projects/brca2/metadata/brca2.metadata.tsv'

###
### 2. annotation
###

#! using biomart: no transcript missing but I end up with 40,320 genes
working_atributes = c('ensembl_transcript_id', 'ensembl_gene_id', 'gene_biotype', 'description')
ensembl96 = useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", version=96)
table96 = getBM(attributes=working_atributes, mart=ensembl96)
dim(table96)
View(table96)

###
### 3. read files
###
metadata = read.table(metadata_file, header=TRUE, sep="\t")
local_samples = metadata$folder
files = file.path(kallisto_dir, local_samples, "abundance.h5")
labels = sapply(strsplit(files, split='/',fixed=TRUE), function(x) (x[9]))
print(files)
print(labels)

txi = tximport(files, type="kallisto", tx2gene=table96, ignoreTxVersion=TRUE)

###
### 4. find abundance
###
tpm = txi$abundance
colnames(tpm) = labels
dim(tpm)
View(tpm)

# 7. store
store = paste(results_dir, 'DESeq2_TPM_values.tsv', sep='')
write.table(tpm, file=store, quote=FALSE, sep='\t', col.names=NA)

