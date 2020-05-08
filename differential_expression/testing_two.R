

# 0. user defined variables
register(MulticoreParam(8))
setwd("~/scratch/")
kallisto_dir = "/Volumes/sand/vigur/data/kallisto_shared_folders"
metadata_file = '/Volumes/sand/vigur/data/metadata/vigur_metadata_experiment_both.tsv'
results_dir = '/Volumes/sand/vigur/results/deseq2/'

###
### 1. read data
###

# 1.1. build annotation reference
txdb = EnsDb.Hsapiens.v86
k = keys(txdb, keytype = "TXNAME")
tx2gene = select(txdb, k, "GENEID", "TXNAME")

# 1.2. read metadata
metadata = read.table(metadata_file, header = TRUE)
rownames(metadata) = metadata$condition
metadata


# 1.3. arrange metadata for hypothesis testing
reference_samples = c(1:3)
hypo = list('reference'=reference_samples, 'tag'='both', 'testing'=c(4:6))

tag = hypo$tag
reference = hypo$reference
testing = hypo$testing

print(tag)
print(reference)
print(testing)

### f.1. define hypothesis metadata
selected_samples = c(reference, testing)
hypothesis_metadata = metadata[selected_samples, ] 
print(hypothesis_metadata)

### f.2. define and read kallisto files
files = file.path(kallisto_dir, hypothesis_metadata$sample, "abundance.h5")
names(files) = hypothesis_metadata$sample
print(files)

txi = tximport(files, type="kallisto", tx2gene=tx2gene, ignoreAfterBar=TRUE, ignoreTxVersion=TRUE)

### f.3. import into DESeq object
dds = DESeqDataSetFromTximport(txi, colData=hypothesis_metadata, design=~ time)
dds$time <- relevel(dds$time, ref = "zero")

### f.4. filtering
threshold = 10
keep = rowSums(counts(dds)) >= threshold  
dds = dds[keep, ]
print(paste('dimensions of filtered genes', dim(dds)[1], dim(dds)[2]))

### f.5 analysis
dds = DESeq(dds, parallel=TRUE)

# 6.2. filter, annotate, format and store

res = results(dds, lfcThreshold=1, parallel=TRUE)
filt1 = res[which(res$pvalue < 0.05), ]
summary(filt1)
filt2 = filt1[which(filt1$padj < 0.1), ]
write.table(filt2, file='testing_c.tsv', quote=FALSE, sep='\t')
