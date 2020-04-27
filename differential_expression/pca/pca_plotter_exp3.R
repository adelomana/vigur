#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("tximport")
#BiocManager::install("EnsDb.Hsapiens.v86")
#BiocManager::install("apeglm")
#BiocManager::install("vsn")
#BiocManager::install("hexbin")
#BiocManager::install("BiocParallel")

library(DESeq2)                
library(tximport)             
library(EnsDb.Hsapiens.v86)
library(apeglm)
library(vsn)
library(hexbin)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(BiocParallel)
library(ggplot2)

# 0. user defined variables
register(MulticoreParam(8))
setwd("~/scratch/")
kallisto_dir = "/Volumes/sand/vigur/data/kallisto_shared_folders"
metadata_file = '/Volumes/sand/vigur/data/metadata/vigur_metadata_experiment3.tsv'

###
### 1. read data
###
# 1.1. build annotation reference
txdb = EnsDb.Hsapiens.v86
k = keys(txdb, keytype = "TXNAME")
tx2gene = select(txdb, k, "GENEID", "TXNAME")
tx2gene

# 1.2. read metadata
endo_metadata = read.table(metadata_file, header = TRUE)
rownames(endo_metadata) <- endo_metadata$condition
endo_metadata

# 1.3. define and read kallisto files
files = file.path(kallisto_dir, endo_metadata$sample, "abundance.h5")
names(files) = endo_metadata$sample
files

txi = tximport(files, type="kallisto", tx2gene=tx2gene, ignoreAfterBar=TRUE, ignoreTxVersion=TRUE)
head(txi$abundance, 5)
head(txi$counts, 5)
dim(txi$counts)

# 1.4. import into DESeq object
dds = DESeqDataSetFromTximport(txi, colData=endo_metadata, design=~ treatment)
dds$treatment <- relevel(dds$treatment, ref = "zero")

# 1.5. filtering
threshold = 100   ### dim: 61,881 --> 20,438
keep = (rowSums(counts(dds))/dim(dds)[2]) >= threshold  
dds = dds[keep,]
dds

###
### 2. analysis
###
dds = DESeq(dds, parallel=TRUE)

###
### 3. sample comparison
###

# 3.1. transformations
vsd = vst(dds, blind=FALSE)

# 3.2. PCA
pcaData = plotPCA(vsd, intgroup=c("treatment", "time"), returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=time)) + 
  scale_color_manual(breaks=c("zero", "half", 'five', 'fifty'), values = c("black", "#59a14f", "#edc948", "#e15759")) + 
  scale_shape_manual(breaks=c("zero", "four", 'twentyfour'), values=c(8, 15, 17)) + 
  geom_point(size=5, alpha=2/3) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()

ggsave('figure_experiment_3.pdf', width=8, height=5, scale=0.8)
