#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("tximport")
#BiocManager::install("EnsDb.Hsapiens.v86")
#BiocManager::install("apeglm")
#BiocManager::install("vsn")
#BiocManager::install("hexbin")
#BiocManager::install("BiocParallel")

# info--
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html bottom

# multiple variables?

library(DESeq2)                
library(tximport)              # to read h5 files faster
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
#metadata_file = '/Volumes/sand/vigur/data/metadata2.txt'
metadata_file = '/Volumes/sand/vigur/data/metadata/vigur_metadata_experiment2.tsv'

# 1. build annotation reference
txdb = EnsDb.Hsapiens.v86
k = keys(txdb, keytype = "TXNAME")
tx2gene = select(txdb, k, "GENEID", "TXNAME")
tx2gene

# 2. read metadata
endo_metadata = read.table(metadata_file, header = TRUE)
rownames(endo_metadata) <- endo_metadata$condition
endo_metadata

# 2. define kallisto files
files = file.path(kallisto_dir, endo_metadata$sample, "abundance.h5")
names(files) = endo_metadata$sample
files

# 3. read kallisto files 
txi = tximport(files, type="kallisto", tx2gene=tx2gene, ignoreAfterBar=TRUE, ignoreTxVersion=TRUE)
head(txi$abundance, 5)
head(txi$counts, 5)

# 1.4. import into DESeq object
dds = DESeqDataSetFromTximport(txi, colData=endo_metadata, design=~ treatment)
#dds = DESeqDataSetFromTximport(txi, colData=endo_metadata, design=~ treatment + time + treatment:time)
dds$treatment <- relevel(dds$treatment, ref = "zero")
#dds$time <- relevel(dds$time, ref = "zero")

# 1.5. filtering
threshold = 1     ### dim: 38640 27 | 349, 198 and 56
threshold = 10    ### dim: 28811 27 | 317, 201 and 55
threshold = 100   ### dim: 17310 27 | 199, 599 and 64
keep = (rowSums(counts(dds))/dim(dds)[2]) >= threshold  
dds = dds[keep,]
dds

# 2. analysis
# 2.1. DESeq
dds = DESeq(dds, parallel=TRUE)
#dds = DESeq(dds, test = "LRT", parallel = TRUE)
res = results(dds, parallel=TRUE)
res

resultsNames(dds)
resA = results(dds, name="treatment_fifty_vs_zero", parallel=TRUE)
resB = results(dds, name="treatment_five_vs_zero", parallel=TRUE)
resC = results(dds, name="treatment_half_vs_zero", parallel=TRUE)

resOrderedA = resA[order(resA$pvalue),]
resOrderedB = resB[order(resB$pvalue),]
resOrderedC = resC[order(resC$pvalue),]

# consider padj < 0.1 and p < 0.05
sum(resA$padj < 0.05, na.rm=TRUE)
sum(resB$padj < 0.05, na.rm=TRUE)
sum(resC$padj < 0.05, na.rm=TRUE)

# 2.2. log fold shrinkage
resLFCA = lfcShrink(dds, coef="treatment_fifty_vs_zero", type="apeglm", parallel=TRUE)
resLFCB = lfcShrink(dds, coef="treatment_five_vs_zero", type="apeglm", parallel=TRUE)
resLFCC = lfcShrink(dds, coef="treatment_half_vs_zero", type="apeglm", parallel=TRUE)
sum(resLFCA$padj < 0.05, na.rm=TRUE)
sum(resLFCB$padj < 0.05, na.rm=TRUE)
sum(resLFCC$padj < 0.05, na.rm=TRUE)

# 3. figures exploring results
# 3.1 MA plots
plotMA(res, ylim=c(-2, 2))

# 3.1 expression of the top differentially expressed gene
plotCounts(dds, gene=which.min(res$padj), intgroup="treatment")

# 4. write results
write.csv(as.data.frame(resOrderedA), file="/Volumes/sand/vigur/results/deseq2/exp2_resOrderedA.csv")
write.csv(as.data.frame(resOrderedB), file="/Volumes/sand/vigur/results/deseq2/exp2_resOrderedB.csv")
write.csv(as.data.frame(resOrderedC), file="/Volumes/sand/vigur/results/deseq2/exp2_resOrderedC.csv")

###
### 5. sample comparison
###

# 5.1. transformations
vsd = vst(dds, blind=FALSE)
rld = rlog(dds, blind=FALSE)
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

# 5.2. heatmap
sampleDists = dist(t(assay(vsd)))
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = vsd$colnames
colnames(sampleDistMatrix) = NULL
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=viridis::viridis(200))

# 5.3. PCA
pcaData = plotPCA(vsd, intgroup=c("treatment", "time"), returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=time)) + 
  scale_color_manual(breaks=c("zero", "half", 'five', 'fifty'), values = c("black", "#59a14f", "#edc948", "#e15759")) + 
  scale_shape_manual(breaks=c("zero", "four", 'twentyfour'), values=c(8, 15, 17)) + 
  geom_point(size=6, alpha=2/3) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()