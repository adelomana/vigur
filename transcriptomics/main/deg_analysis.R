# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("XML")
# BiocManager::install("openssl")
# BiocManager::install("GenomicFeatures")
# BiocManager::install("DESeq2")
# BiocManager::install("tximport")
# BiocManager::install("EnsDb.Hsapiens.v86")
# BiocManager::install("BiocParallel")
# BiocManager::install("rlist")
# BiocManager::install("stringr")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("tictoc")
# BiocManager::install("rhdf5")

library(DESeq2)                
library(tximport)         
library(EnsDb.Hsapiens.v86)
library(BiocParallel)
library(rlist)
library(stringr)
library(org.Hs.eg.db)
library(tictoc)
library(rhdf5)

tic()

# 0. user defined variables
register(MulticoreParam(8))
setwd("~/scratch/")
kallisto_dir = "/home/adrian/projects/vigur/data/transcriptomics/kallisto_shared_folders"
metadata_file = '/home/adrian/projects/vigur/data/transcriptomics/metadata/vigur_metadata_experiment3.tsv'
results_dir = '/home/adrian/projects/vigur/results/transcriptomics/deseq2/'

###
### 1. read data
###

# 1.1. build annotation reference
k = keys(EnsDb.Hsapiens.v86, keytype = "TXNAME")
tx2gene = select(EnsDb.Hsapiens.v86, k, "GENEID", "TXNAME")

# 1.2. read metadata
metadata = read.table(metadata_file, header = TRUE)
rownames(metadata) = metadata$condition
metadata

# 2. store abundance for all samples
files = file.path(kallisto_dir, metadata$sample, "abundance.h5")
txi = tximport(files, type="kallisto", tx2gene=tx2gene, ignoreAfterBar=TRUE, ignoreTxVersion=TRUE)

tpm = txi$abundance
colnames(tpm) = metadata$sample
dim(tpm)

ensembl_ids = rownames(tpm)
annotation = select(EnsDb.Hsapiens.v86, ensembl_ids, c("GENEBIOTYPE", "GENENAME"), "GENEID")

store = paste(results_dir, 'DESeq2_TPM_values.tsv', sep='')
write.table(tpm, file=store, quote=FALSE, sep='\t', col.names=NA)

store = paste(results_dir, 'annotation.tsv', sep='')
write.table(annotation, file=store, quote=FALSE, sep='\t', col.names=NA)

# 3. arrange metadata for hypothesis testing
hypos = list()

### experiment 3

# concentration zero
reference_samples = c(1:3) # time zero, concentration zero
hypo = list('reference'=reference_samples, 'tag'='experiment_three_concentration_zero_time_four_vs_control', 'testing'=c(4:6), 'influence'='time')
hypos = list.append(hypos, hypo)
hypo = list('reference'=reference_samples, 'tag'='experiment_three_concentration_zero_time_twentyfour_vs_control', 'testing'=c(14:16), 'influence'='time')
hypos = list.append(hypos, hypo)

# concentration half
reference_samples = c(4:6) # time four, concentration zero
hypo = list('reference'=reference_samples, 'tag'='experiment_three_concentration_half_time_four_vs_control', 'testing'=c(7:8), 'influence'='concentration')
hypos = list.append(hypos, hypo)
reference_samples = c(14:16) # time twentyfour, concentration zero
hypo = list('reference'=reference_samples, 'tag'='experiment_three_concentration_half_time_twentyfour_vs_control', 'testing'=c(17:19), 'influence'='concentration')
hypos = list.append(hypos, hypo)

# concentration five
reference_samples = c(4:6) # time four, concentration zero
hypo = list('reference'=reference_samples, 'tag'='experiment_three_concentration_five_time_four_vs_control', 'testing'=c(9:11), 'influence'='concentration')
hypos = list.append(hypos, hypo)
reference_samples = c(14:16) # time twentyfour, concentration zero
hypo = list('reference'=reference_samples, 'tag'='experiment_three_concentration_five_time_twentyfour_vs_control', 'testing'=c(20:22), 'influence'='concentration')
hypos = list.append(hypos, hypo)

# concentration fifty
reference_samples = c(4:6) # time four, concentration zero
hypo = list('reference'=reference_samples, 'tag'='experiment_three_concentration_fifty_time_four_vs_control', 'testing'=c(12:13), 'influence'='concentration')
hypos = list.append(hypos, hypo)
reference_samples = c(14:16) # time twentyfour, concentration zero
hypo = list('reference'=reference_samples, 'tag'='experiment_three_concentration_fifty_time_twentyfour_vs_control', 'testing'=c(23:25), 'influence'='concentration')
hypos = list.append(hypos, hypo)

# 4. define analysis function
compare = function(hypo) {
  
  tag = hypo$tag
  reference = hypo$reference
  testing = hypo$testing
  
  print('---------- working on a new hypothesis')
  print(tag)
  print('reference')
  print(reference)
  print('testing samples')
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
  if (hypo$influence == 'time'){
    print('working with time as influence...')
    dds = DESeqDataSetFromTximport(txi, colData=hypothesis_metadata, design=~time)
    dds$time <- relevel(dds$time, ref="zero")
  }
  if (hypo$influence == 'concentration'){
    print('working with concentration as influence...')
    dds = DESeqDataSetFromTximport(txi, colData=hypothesis_metadata, design=~treatment)
    dds$treatment <- relevel(dds$treatment, ref="zero")
  }
  
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
  print(paste('DEGs found', dim(filt2)[1], sep=' '))
  write.table(filt2, file=paste(results_dir, 'unformatted/unformatted_results.', tag, '.tsv', sep=''), quote=FALSE, sep='\t')
  
  print('annotate')
  df = as.data.frame(filt2)
  df['common'] = rownames(df)
  selected = rownames(df)
  info = select(EnsDb.Hsapiens.v86, selected, c("GENEBIOTYPE", "GENENAME"), "GENEID")
  info['common'] = info$GENEID
  
  descriptions = tryCatch({
    descriptions = select(org.Hs.eg.db, keys=selected, columns=c("GENENAME"), keytype="ENSEMBL")
  }, error = function(e) {
    print('Warning: no description found for ENSEMBL IDs')
    descriptions = data.frame('ENSEMBL'=selected, 'GENENAME'=rep('Not found', each=length(selected)))
  })
  
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
  store = paste(results_dir, tag, '_up', '.tsv', sep='')
  print(store)
  write.table(sorted_up, file=store, quote=FALSE, sep='\t', row.names=FALSE)
  store = paste(results_dir, tag, '_down', '.tsv', sep='')
  print(store)
  write.table(sorted_down, file=store, quote=FALSE, sep='\t', row.names=FALSE)
  
}

# 5. iterate function
for (hypo in hypos) {
  tempo = compare(hypo)
}

toc()
