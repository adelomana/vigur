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
results_dir = '/home/adrian/projects/vigur/results/sequencing/tpm/'
kallisto_dir = "/home/adrian/projects/vigur/results/sequencing/kallisto/kallisto.100"

###
### 2. annotation
###

#! using file provided by kallisto, I get transcripts missing from tx2gene: 17686 and then a dim(tpm) - 35,606
#tx2gene_map = read.csv(annotation_file, sep='\t', header=FALSE, col.names=c('TXID', 'GENEID', 'gene_name'))
#tx2gene_map

#! using v86 package
#! transcripts missing from tx2gene: 11043
#! i end up with a tpm matrix of 39,415
# library(EnsDb.Hsapiens.v86)
# k = keys(EnsDb.Hsapiens.v86, keytype = "TXNAME")
# tx2gene86 = select(EnsDb.Hsapiens.v86, k, "GENEID", "TXNAME")
# view(tx2gene86)

#! using biomart: no transcript missing but I end up with 40,320 genes
working_atributes = c('ensembl_transcript_id', 'ensembl_gene_id', 'gene_biotype', 'description')
ensembl96 = useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", version=96)
table96 = getBM(attributes=working_atributes, mart=ensembl96)
dim(table96)
View(table96)
ordered_table96 = table96[, c(2, 1, 3, 4)]
dim(ordered_table96)
View(ordered_table96)

###
### 3. read files
###
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
dirnames

files = file.path(dirnames, 'abundance.h5')
print(files)
labels = sapply(strsplit(files, split='/',fixed=TRUE), function(x) (x[10]))

txi = tximport(files, type="kallisto", tx2gene=ordered_table96, ignoreTxVersion=TRUE)

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

store = paste(results_dir, 'annotation.tsv', sep='')
write.table(ordered_table96, file=store, quote=FALSE, sep='\t', col.names=NA)

