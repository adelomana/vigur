if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "tximport", "EnsDb.Hsapiens.v86"))

library('DESeq2')
library("tximport")
library("EnsDb.Hsapiens.v86")

# 0. user defined variables
setwd("~/scratch/")
kallisto_dir = "/Volumes/sand/vigur/data/kallisto_shared_folders"
metadata_file = '/Volumes/sand/vigur/data/metadata.txt'



# 1. read metadata
endo_metadata = read.table(metadata_file, header = TRUE)
endo_metadata



##### annotation
txdb <- EnsDb.Hsapiens.v86
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
tx2geje

# 2. read kallisto files
files = file.path(kallisto_dir, endo_metadata$run, "abundance.tsv")
files

# 2. 
txi = tximport(files, type="salmon", tx2gene=tx2gene)
ddsTxi = DESeqDataSetFromTximport(txi, colData=samples, design = ~ condition)