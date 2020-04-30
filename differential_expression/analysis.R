#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")

library(DESeq2)                
library(tximport)         
library(EnsDb.Hsapiens.v86)
library(BiocParallel)
library(rlist)
library(stringr)
library(org.Hs.eg.db)

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

# 2. store abundance for all samples
#tpm <- txi$abundance
#store = paste(results_dir, 'DESeq2_TPM_values.tsv', '.tsv', sep='')
#write.csv(tpm, file=store, quote=FALSE, sep='\t')

# 3. arrange metadata for hypothesis testing
hypos = list()

### experiment 2
reference_samples = c(1:3)

# concentration zero
hypo = list('reference'=reference_samples, 'tag'='experiment_two_concentration_zero_time_four_vs_zero', 'testing'=c(16:18))
hypos = list.append(hypos, hypo)
hypo = list('reference'=reference_samples, 'tag'='experiment_two_concentration_zero_time_twentyfour_vs_zero', 'testing'=c(4:6))
hypos = list.append(hypos, hypo)
#c(4:6)
#c(16:18)

# concentration half
hypo = list('reference'=reference_samples, 'tag'='experiment_two_concentration_half', 'testing'=c(c(7:9), c(19:21)))
#hypos = list.append(hypos, hypo)

# concentration five
hypo = list('reference'=reference_samples, 'tag'='experiment_two_concentration_five', 'testing'=c(c(10:12), c(22:24)))
#hypos = list.append(hypos, hypo)

# concentration fifty
hypo = list('reference'=reference_samples, 'tag'='experiment_two_concentration_fifty', 'testing'=c(c(13:15), c(25:27)))
#hypos = list.append(hypos, hypo)

### experiment 3
reference_samples = c(28:30)

# concentration zero
hypo = list('reference'=reference_samples, 'tag'='experiment_three_concentration_zero', 'testing'=c(c(31:33), c(41:43)))
#hypos = list.append(hypos, hypo)

# concentration half
hypo = list('reference'=reference_samples, 'tag'='experiment_three_concentration_half', 'testing'=c(c(34:35), c(44:46)))
#hypos = list.append(hypos, hypo)

# concentration five
hypo = list('reference'=reference_samples, 'tag'='experiment_three_concentration_five', 'testing'=c(c(36:38), c(47:49)))
#hypos = list.append(hypos, hypo)

# concentration fifty
hypo = list('reference'=reference_samples, 'tag'='experiment_three_concentration_fifty', 'testing'=c(c(39:40), c(50:52)))
#hypos = list.append(hypos, hypo)

# 4. define analysis function
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
  
  ### f.3. import into DESeq object
  dds = DESeqDataSetFromTximport(txi, colData=hypothesis_metadata, design=~time)
  dds$time <- relevel(dds$time, ref = "zero")
  
  ### f.4. filtering
  threshold = 10
  keep = rowSums(counts(dds)) >= threshold  
  dds = dds[keep, ]
  print(paste('dimensions of filtered genes', dim(dds)[1], dim(dds)[2]))
  
  ### f.5 analysis
  print('analysis')
  dds = DESeq(dds, parallel=TRUE)

  # f.6. filter, annotate, format and store
  print('filter')
  res = results(dds, lfcThreshold=1, parallel=TRUE)
  filt1 = res[which(res$pvalue < 0.05), ]
  filt2 = filt1[which(filt1$padj < 0.1), ]
  write.table(filt2, file='testing.tsv', quote=FALSE, sep='\t')
  
  print('annotate')
  df = as.data.frame(filt2)
  df['common'] = rownames(df)
  selected = rownames(df)
  info = as.data.frame(select(txdb, selected, c("GENEBIOTYPE", "GENENAME"), "GENEID"))
  info['common'] = info$GENEID
  print('before')
  descriptions = as.data.frame(select(org.Hs.eg.db, keys=selected, columns=c("GENENAME"), keytype="ENSEMBL"))
  print('after')
  names(descriptions)[names(descriptions) == "GENENAME"] <- "DESCRIPTION" # arrow is needed here!
  descriptions['common'] = descriptions$ENSEMBL
  dh = merge(df, info, by='common')
  di = merge(dh, descriptions, by='common')
  
  print('format')
  formatted = di[ , c(8,10,9,12,2,3,6,7)]
  up = formatted[formatted$log2FoldChange > 0, ]
  down = formatted[formatted$log2FoldChange < 0, ]
  sorted_up = up[order(up$log2FoldChange, decreasing=TRUE), ]
  sorted_down = down[order(down$log2FoldChange), ]
  
  print('store')
  store = paste(results_dir, tag, '_up_', '.tsv', sep='')
  write.table(sorted_up, file=store, quote=FALSE, sep='\t', row.names=FALSE)
  store = paste(results_dir, tag, '_down_', '.tsv', sep='')
  write.table(sorted_down, file=store, quote=FALSE, sep='\t', row.names=FALSE)
  
  print('---')
}

# 5. iterate function
for (hypo in hypos) {
  tempo = compare(hypo)
}
