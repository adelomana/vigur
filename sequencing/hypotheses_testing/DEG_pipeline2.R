###
### This script takes all six samples, mising biological and technical replicates and finds DEGs over all treatment conditions for both 072 and 073.
### Sample 073_39 (mixed cat + ilo @ 4h, rep 2) is explicitly excluded because the sequencing failed.
###

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("tictoc")

###
### 1. load libraries
###
library(rlist) # necessary for list appending
library(biomaRt)
library(tximport)
library(DESeq2)

# performance
library(BiocParallel)
library(tictoc)

###
### 2. user defined variables
###
register(MulticoreParam(20))

setwd("~/scratch/")
kallisto_dir = "/home/adrian/projects/vigur/results/sequencing/kallisto/kallisto.100"
metadata_file = '/home/adrian/projects/vigur/data/sequencing/metadata/metadata.tsv'
results_dir = '/home/adrian/projects/vigur/results/sequencing/DEGs/'

###
### 3. read gene annotation
###
#! using biomart, fuck yeah
working_atributes = c('ensembl_transcript_id', 'ensembl_gene_id', 'gene_biotype', 'description', 'hgnc_symbol')
ensembl96 = useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", version=96)
#! listAttributes(mart=ensembl96)
table96 = getBM(attributes=working_atributes, mart=ensembl96)
table96['common'] = table96$ensembl_gene_id
#dim(table96)
#View(table96)
#! the order should be first transcript as first column and gene as second column. The rest of columns are not used.

###
### 3. read metadata
###
metadata = read.table(metadata_file, header=TRUE, sep="\t")
rownames(metadata) = metadata$sampleID # this line is necessary for later local selection of metadata during hypothesis testing

###
### 4. define hypotheses
###
hypos = list() # beware this is a list, not a vector as c()

###
### 4.1. define hypotheses for 072
###
run = 72
timepoints = c('four', 'twentyfour')
treatments = c('five_epi', 'five_nor', 'mix', 'TNFa')

### define first the time contrasts without treatments
influence = 'time'

reference_samples = metadata[metadata$run == run & metadata$treatment == 'zero' & metadata$time == 'zero', ]$sampleID
testing_samples = metadata[metadata$run == run & metadata$treatment == 'zero' & metadata$time == 'four', ]$sampleID

print('testing')
print(metadata[which(metadata$sampleID %in% testing_samples), ])
print('reference')
print(metadata[which(metadata$sampleID %in% reference_samples), ])

the_tag = paste('run', run, influence, 'four_vs_zero', sep='_')
hypo = list('reference'=reference_samples, 'testing'=testing_samples, 'influence'=influence, 'tag'=the_tag)
print(hypo)
hypos = list.append(hypos, hypo)
print('----------------------------------------------------------')

reference_samples = metadata[metadata$run == run & metadata$treatment == 'zero' & metadata$time == 'zero', ]$sampleID
testing_samples = metadata[metadata$run == run & metadata$treatment == 'zero' & metadata$time == 'twentyfour', ]$sampleID

print('testing')
print(metadata[which(metadata$sampleID %in% testing_samples), ])
print('reference')
print(metadata[which(metadata$sampleID %in% reference_samples), ])

the_tag = paste('run', run, influence, 'twentyfour_vs_zero', sep='_')
hypo = list('reference'=reference_samples, 'testing'=testing_samples, 'influence'=influence, 'tag'=the_tag)
print(hypo)
hypos = list.append(hypos, hypo)
print('----------------------------------------------------------')

##############################################################

### working with treatment contrasts
influence = 'treatment'

for (time in timepoints){
  for (treatment in treatments){
    reference_samples = metadata[metadata$run == run & metadata$treatment == 'zero' & metadata$time == time, ]$sampleID
    testing_samples = metadata[metadata$run == run & metadata$treatment == treatment & metadata$time == time, ]$sampleID
    
    print('testing')
    print(metadata[which(metadata$sampleID %in% testing_samples), ])
    print('reference')
    print(metadata[which(metadata$sampleID %in% reference_samples), ])
    
    the_tag = paste('run', run, 'treatment', treatment, 'time', time, sep='_')
    hypo = list('reference'=reference_samples, 'testing'=testing_samples, 'influence'=influence, 'tag'=the_tag)
    print(hypo)
    hypos = list.append(hypos, hypo)
    print('----------------------------------------------------------')
    
  }
}

###
### 4.2. define hypotheses for 073
###
run = 73
timepoints = c('four')
treatments = c('mix', 'ilo_only', 'mix_plus_ilo')

### define first the time contrasts without treatments
influence = 'time'

reference_samples = metadata[metadata$run == run & metadata$treatment == 'zero' & metadata$time == 'zero', ]$sampleID
testing_samples = metadata[metadata$run == run & metadata$treatment == 'zero' & metadata$time == 'four', ]$sampleID

print('testing')
print(metadata[which(metadata$sampleID %in% testing_samples), ])
print('reference')
print(metadata[which(metadata$sampleID %in% reference_samples), ])

the_tag = paste('run', run, influence, 'four_vs_zero', sep='_')
hypo = list('reference'=reference_samples, 'testing'=testing_samples, 'influence'=influence, 'tag'=the_tag)
print(hypo)
hypos = list.append(hypos, hypo)
print('----------------------------------------------------------')

##############################################################

### working with treatment contrasts
influence = 'treatment'

for (time in timepoints){
  for (treatment in treatments){
    reference_samples = metadata[metadata$run == run & metadata$treatment == 'zero' & metadata$time == time, ]$sampleID
    testing_samples = metadata[metadata$run == run & metadata$treatment == treatment & metadata$time == time, ]$sampleID
    
    # making sure that xx is not used
    testing_samples = testing_samples[testing_samples!= '073_39']
    
    print('testing')
    print(metadata[which(metadata$sampleID %in% testing_samples), ])
    print('reference')
    print(metadata[which(metadata$sampleID %in% reference_samples), ])
    
    the_tag = paste('run', run, 'treatment', treatment, 'time', time, sep='_')
    hypo = list('reference'=reference_samples, 'testing'=testing_samples, 'influence'=influence, 'tag'=the_tag)
    print(hypo)
    hypos = list.append(hypos, hypo)
    print('----------------------------------------------------------')
    
  }
}

###
### 5. run the hypotheses
###
#hypos = c()
tic()
for (hypo in hypos) {
  
  tag = hypo$tag
  influence = hypo$influence
  reference = hypo$reference
  testing = hypo$testing
  local_samples = c(reference, testing)
  local_metadata = metadata[local_samples, ]
  
  print('------------ working on a new hypothesis')
  print(tag)
  print('reference')
  print(reference)
  print('testing samples')
  print(testing)
  print('influence')
  print(influence)
  print('local samples')
  print(local_samples)
  print('local metadata')
  print(local_metadata)
  
  # 5.1. define and read files
  files = file.path(kallisto_dir, local_samples, "abundance.h5")
  print(files)
  txi = tximport(files, type="kallisto", tx2gene=table96, ignoreTxVersion=TRUE)
  
  # 5.2. define DESeq2 object
  if (influence == 'time'){
    print('working with time as influence...')
    dds = DESeqDataSetFromTximport(txi, colData=local_metadata, design=~time) 
    dds$time <- relevel(dds$time, ref="zero")
  }
  if (influence == 'treatment'){
    print('working with concentration as influence...')
    dds = DESeqDataSetFromTximport(txi, colData=local_metadata, design=~treatment) 
    dds$treatment <- relevel(dds$treatment, ref="zero")
  }
  
  ### 5.3. filtering
  threshold = 5
  keep = rowMaxs(counts(dds)) >= threshold
  dds = dds[keep, ]
  print(paste('dimensions of filtered genes', dim(dds)[1], dim(dds)[2]))
  
  ### 5.4. analysis
  print('analysis')
  dds = DESeq(dds, parallel=TRUE)
  
  ### 5.5. filter, annotate, format and store
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
  #View(annotation_table_unique)
  
  dn = merge(df, annotation_table_unique, by='common')
  #View(dn)
  
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
  #View(sorted_up)
  
  print('store')
  store = paste(results_dir, tag, '_up', '.tsv', sep='')
  print(store)
  write.table(sorted_up, file=store, quote=FALSE, sep='\t', row.names=FALSE)
  store = paste(results_dir, tag, '_down', '.tsv', sep='')
  print(store)
  write.table(sorted_down, file=store, quote=FALSE, sep='\t', row.names=FALSE)
  
}
toc()










  