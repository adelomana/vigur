#!/usr/bin/env Rscript

##
### install IPO
###

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# # IPO dependencies. They may vary depending on your machine
# BiocManager::install("foreign")
# BiocManager::install("spatial")
# BiocManager::install("igraph")
# 
# BiocManager::install('XML')
# BiocManager::install('mzID') 
# BiocManager::install('MSnbase')
# BiocManager::install('mzR')
# BiocManager::install('SummarizedExperiment')
# 
# BiocManager::install("CAMERA")
# BiocManager::install("IPO")

###
### load libraries
###

library(IPO)
library(BiocParallel)

###
### 0. user-defined variables
###

# paths
data_dir = '/home/adrian/projects/vigur/data/metabolomics/pools/'
results_dir = '/home/adrian/projects/vigur/results/metabolomics/pools/'

# variables
iters = 1
nThreads = 20

### 1. read info
cases = list.files(data_dir)
cases = c("negative", "positive")

###
### 2. analysis
###
for (case_index in 1:length(cases)) {
  label = cases[case_index]
  results_subdir = paste(results_dir, label, sep='')
  dir.create(results_subdir)

  subdir = paste(data_dir, label, sep='')
  datafiles = list.files(subdir, recursive=TRUE, full.names=TRUE) 
  print(datafiles)
  
  for (iter in 1:iters) {
    
    iter_subdir = paste(results_subdir, iter, sep='/')
    dir.create(iter_subdir)
    setwd(iter_subdir)
    
    # explore parameters
    ppParameters = getDefaultXcmsSetStartingParams('centWave')
    
    # freyrs boundaries
    #ppParameters$min_peakwidth = c(2, 7)
    #ppParameters$max_peakwidth = c(10, 20)
    #ppParameters$ppm = c(5, 50)
    #ppParameters$mzdiff = c(-0.05, 0.010)
    
    resultParameters = optimizeXcmsSet(files=datafiles, params=ppParameters, plot=TRUE, BPPARAM=MulticoreParam(workers=nThreads))
      
    # write results
    results_file = paste(label, iter, '.tsv', sep='')
    df = as.data.frame(resultParameters$best_settings$parameters)
    write.table(df, file=results_file, sep='\t', quote=FALSE)
  }
}
