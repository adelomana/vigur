###
### load libraries
###

library(IPO)
library(BiocParallel)

###
### 0. user-defined variables
###

# paths
data_dir = '/home/adrian/projects/vigur/data/metabolomics/pools/mini/'
results_dir = '/home/adrian/projects/vigur/results/metabolomics/pools/'

# variables
nThreads = 2

###
### 2. analysis
###
label = 'mini'
results_subdir = paste(results_dir, label, sep='')
dir.create(results_subdir)

datafiles = list.files(data_dir, recursive=TRUE, full.names=TRUE) 
print(datafiles)
  
setwd(results_subdir)
    
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