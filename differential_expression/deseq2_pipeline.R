#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("DESeq2", "tximport", "EnsDb.Hsapiens.v86"))

library('DESeq2')
library("tximport")
library("EnsDb.Hsapiens.v86")

# 0. user defined variables
setwd("~/scratch/")
kallisto_dir = "/Volumes/sand/vigur/data/kallisto_shared_folders"
metadata_file = '/Volumes/sand/vigur/data/metadata.txt'

# 1. build annotation reference
txdb = EnsDb.Hsapiens.v86
k = keys(txdb, keytype = "TXNAME")
tx2gene = select(txdb, k, "GENEID", "TXNAME")
tx2geje

# 2. read metadata
endo_metadata = read.table(metadata_file, header = TRUE)
endo_metadata$condition <- factor(rep(c("A","B"),each=3))
rownames(endo_metadata) <- endo_metadata$run
endo_metadata

# 2. define kallisto files
files = file.path(kallisto_dir, endo_metadata$run, "abundance.h5")
files
names(files) <- paste0("sample", 1:6)

# 3. read kallisto files 
txi = tximport(files, type="kallisto", tx2gene=tx2gene, ignoreAfterBar=TRUE, ignoreTxVersion=TRUE)
head(txi$counts)

# analysis
ddsTxi = DESeqDataSetFromTximport(txi, colData=metadata, design=~condition)
