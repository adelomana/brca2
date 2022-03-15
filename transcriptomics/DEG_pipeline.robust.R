###
### This script compares two groups: +/- and -/-
###

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("tictoc")

###
### 1. load libraries
###
library(biomaRt)
library(tximport)
library(DESeq2)

# performance
library(BiocParallel)

###
### 2. user defined variables
###
register(MulticoreParam(20))

setwd("~/scratch/")
kallisto_dir = "/home/adrian/projects/brca2/results/kallisto/kallisto.100"
metadata_file = '/home/adrian/projects/brca2/metadata/metadata_verbose.tsv'
results_dir = '/home/adrian/projects/brca2/results/DESeq2/'

###
### 3. read gene annotation
###
#! using biomart, fuck yeah
working_atributes = c('ensembl_transcript_id', 'ensembl_gene_id', 'gene_biotype', 'description', 'hgnc_symbol')
ensembl96 = useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", version=96)
#! listAttributes(mart=ensembl96)
table96 = getBM(attributes=working_atributes, mart=ensembl96)
table96['common'] = table96$ensembl_gene_id
dim(table96)
View(table96)
#! the order should be first transcript as first column and gene as second column. The rest of columns are not used.

###
### 3. read metadata
###
metadata = read.table(metadata_file, header=TRUE, sep="\t")
metadata = metadata[-1,]
View(metadata)

###
### 4. run contrasts
###

tag = 'robust_comparison'
local_samples = metadata$true.labels

# 4.1. define and read files
files = file.path(kallisto_dir, local_samples, "abundance.h5")
print(files)
txi = tximport(files, type="kallisto", tx2gene=table96, ignoreTxVersion=TRUE)
  
# 4.2. define DESeq2 object
dds = DESeqDataSetFromTximport(txi, colData=metadata, design=~condition) 
dds$condition <- relevel(dds$condition, ref="neg")
print(paste('dimensions of unfiltered genes', dim(dds)[1], dim(dds)[2]))
 
# 4.3. filtering
threshold = 5
keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]
print(paste('dimensions of filtered genes', dim(dds)[1], dim(dds)[2]))
  
# 4.4. analysis
print('analysis')
dds = DESeq(dds, parallel=TRUE)
  
# 5.5. filter, annotate, format and store
print('filter')
res = results(dds, lfcThreshold=1, parallel=TRUE) 
filt1 = res[which(res$pvalue < 0.05), ]
filt2 = filt1[which(filt1$padj < 0.1), ]
print(paste('DEGs found', dim(filt2)[1], sep=' '))
write.table(filt2, file=paste(results_dir, 'unformatted/unformatted_results_', tag, '.tsv', sep=''), quote=FALSE, sep='\t')
  
print('annotate')
df = as.data.frame(filt2)
df['common'] = rownames(df)
  
annotation_table = table96[, c(3, 4, 5, 6)]
annotation_table_unique = annotation_table[!duplicated(annotation_table$common), ]
View(annotation_table_unique)
  
dn = merge(df, annotation_table_unique, by='common')
View(dn)
  
# check for not missing DEGs because of annotation
if (dim(df)[1] != dim(dn)[1]){
  print('ERROR: DEG missed on annotation step')
  stop()
}
  
print('format')
up = dn[dn$log2FoldChange > 0, ]
down = dn[dn$log2FoldChange < 0, ]
sorted_up = up[order(up$log2FoldChange, decreasing=TRUE), ]
sorted_down = down[order(down$log2FoldChange), ]
View(sorted_up)
  
print('store')
store = paste(results_dir, tag, '_up', '.tsv', sep='')
print(store)
write.table(sorted_up, file=store, quote=FALSE, sep='\t', row.names=FALSE)
store = paste(results_dir, tag, '_down', '.tsv', sep='')
print(store)
write.table(sorted_down, file=store, quote=FALSE, sep='\t', row.names=FALSE)










  