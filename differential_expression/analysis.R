#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")

library(DESeq2)                
library(tximport)         
library(EnsDb.Hsapiens.v86)
library(BiocParallel)
library(rlist)
library(stringr)

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
hypos = list()

### experiment 2
reference_samples = c(1:3)

# concentration zero
hypo = list('reference'=reference_samples, 'tag'='experiment_two_concentration_zero', 'testing'=c(c(4:6), c(16:18)))
hypos = list.append(hypos, hypo)

# concentration half
hypo = list('reference'=reference_samples, 'tag'='experiment_two_concentration_half', 'testing'=c(c(7:9), c(19:21)))
hypos = list.append(hypos, hypo)

# concentration five
hypo = list('reference'=reference_samples, 'tag'='experiment_two_concentration_five', 'testing'=c(c(10:12), c(22:24)))
hypos = list.append(hypos, hypo)

# concentration fifty
hypo = list('reference'=reference_samples, 'tag'='experiment_two_concentration_fifty', 'testing'=c(c(13:15), c(25:27)))
hypos = list.append(hypos, hypo)

### experiment 3
reference_samples = c(28:30)

# concentration zero
hypo = list('reference'=reference_samples, 'tag'='experiment_three_concentration_zero', 'testing'=c(c(31:33), c(41:43)))
hypos = list.append(hypos, hypo)

# concentration half
hypo = list('reference'=reference_samples, 'tag'='experiment_three_concentration_half', 'testing'=c(c(34:35), c(44:46)))
hypos = list.append(hypos, hypo)

# concentration five
hypo = list('reference'=reference_samples, 'tag'='experiment_three_concentration_five', 'testing'=c(c(36:38), c(47:49)))
hypos = list.append(hypos, hypo)

# concentration fifty
hypo = list('reference'=reference_samples, 'tag'='experiment_three_concentration_fifty', 'testing'=c(c(39:40), c(50:52)))
hypos = list.append(hypos, hypo)
# 2. define analysis function
compare = function(hypo) {
  
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
  print(head(txi$abundance, 3))
  print(head(txi$counts, 3))
  
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
  res = results(dds, parallel=TRUE)
  
  res_four = results(dds, name="time_four_vs_zero", parallel=TRUE)
  res_ordered_four = res_four[order(res_four$padj), ]
  a = sum(res_ordered_four$padj < 0.1, na.rm=TRUE)
  print(paste('T4 DEGs', a))
  
  res_twentyfour = results(dds, name="time_twentyfour_vs_zero", parallel=TRUE)
  res_ordered_twentyfour = res_twentyfour[order(res_twentyfour$padj), ]
  b = sum(res_ordered_twentyfour$padj < 0.1, na.rm=TRUE)
  print(paste('T24 DEGs', b))
  
  ### f.6. store results
  four_file = paste(results_dir, tag, '_', 'time_four_vs_zero', '.csv', sep='')
  twentyfour_file = paste(results_dir, tag, '_', 'time_twentyfour_vs_zero', '.csv', sep='')
  write.csv(as.data.frame(res_ordered_four), file=four_file)
  write.csv(as.data.frame(res_ordered_twentyfour), file=twentyfour_file)
  
  print('---')
}

# 3. iterate function
for (hypo in hypos) {
  tempo = compare(hypo)
}
