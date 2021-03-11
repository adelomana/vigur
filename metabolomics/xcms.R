#!/usr/bin/env Rscript

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("tictoc")

library('IPO')
library('xcms')
library('tictoc')

###
### 0. user-defined variables
###

### 0.1. paths
data_dir = '/home/adrian/projects/vigur/data/metabolomics/sets/'
ipo_optimal_paramters_file = '/home/adrian/projects/vigur/results/metabolomics/pools/ipo_output.RData'
results_dir = '/home/adrian/projects/vigur/results/metabolomics/xcms/'

# 0.2. parameters
nThreads = 20

###
### 1. analysis
###

### 1.1. load IPO parameters
load(ipo_optimal_paramters_file) ### optimal parameters are stored into a variable called "full_set_parameters"
full_set_parameters

###
### 2. xcms iteratations over batches
###

xcms_optimized_parameters = data.frame()

for (batch_index in c(1, 2)){
  message(batch_index)
}

#for (batch_index in 1:nrow(full_set_parameters)){
for (batch_index in c(1, 3)){

  ### 2.1. set batch-specific parameters
  min_pw = full_set_parameters[batch_index, "min_peakwidth"]
  max_pw = full_set_parameters[batch_index, "min_peakwidth"]
  batch_ppm = full_set_parameters[batch_index, "ppm"]
  batch_mzdiff = full_set_parameters[batch_index, "ppm"]
  
  batch_name = rownames(full_set_parameters)[batch_index]
  working_dir = paste(data_dir, batch_name, '_dir', sep='')
  message(working_dir)
  message(batch_name)
  
  if (grepl('positive', batch_name, fixed = TRUE) == TRUE) {batch_polarity = 'positive'} else {batch_polarity = 'negative'} ### basic mode runs at negative polarity
  
  ### 2.2. set the object
  hlutur = xcmsSet(working_dir,  
                   peakwidth=c(min_pw, max_pw),
                   ppm=batch_ppm,
                   mzdiff=batch_mzdiff, 
                   polarity=batch_polarity,
                   BPPARAM=MulticoreParam(workers=nThreads),
                   #
                   prefilter=c(3, 100),
                   method='centWave',
                   snthresh=10, 
                   mzCenterFun="wMean", 
                   integrate=1,
                   noise=0, 
                   fitgauss=FALSE)
  
  ### 2.3. define some parameters
  ud_minfrac = (6-1)/36 # at least 6 samples are biologically similar. 30 samples and 6 pools, 36 in total
  ud_minfrac
  ud_minsamp = 1
  ud_max = 50
  
  retcor_boundaries = getDefaultRetGroupStartingParams()
  message()
  message('original boundaries')
  print(retcor_boundaries)
  
  retcor_boundaries$bw = c(0.1, 10) #  22-38 for HPLC, lower bounds for UPLC, which we are using. Freyrs bounds were (4, 25)
  retcor_boundaries$minfrac = ud_minfrac
  retcor_boundaries$mzwid = c(0.005, 0.2)
  retcor_boundaries$profStep = c(0.05, 1)
  retcor_boundaries$gapInit = c(0.05, 2)
  #retcor_boundaries$gapExtend = c(1, 50)
  message()
  message('modified boundaries')
  print(retcor_boundaries)
  
  optimized_retcor_p = optimizeRetGroup(xset=hlutur, params=retcor_boundaries, nSlaves=nThreads)
  
  print('Optimization results:')
  print(c('bw', retcor_boundaries$bw, optimized_retcor_p$best_settings$bw)) # (2,10) 2   
  print(c('minfrac', retcor_boundaries$minfrac, optimized_retcor_p$best_settings$minfrac)) # (0.05, 0.5) 0.05
  print(c('mzwid', retcor_boundaries$mzwid, optimized_retcor_p$best_settings$mzwid)) # 0.015 0.035 0.015
  print(c('profStep', retcor_boundaries$profStep, optimized_retcor_p$best_settings$profStep)) # 0.5 1 0.5
  print(c('gapInit', retcor_boundaries$gapInit, optimized_retcor_p$best_settings$gapInit)) # 0 0.8 0.8
  print(c('gapExtend', retcor_boundaries$gapExtend, optimized_retcor_p$best_settings$gapExtend)) # 2.1 2.7 2.4
  
  optimized_retcor_p$best_settings$label = batch_name
  optimized_retcor_p$best_settings
  
  message('about to add xcms optimized parameters...')
  xcms_optimized_parameters = rbind(xcms_optimized_parameters, optimized_retcor_p$best_settings)
  
  message('about one')
  grouped = group(hlutur, bw=optimized_retcor_p$best_settings$bw, minfrac=ud_minfrac, mzwid=optimized_retcor_p$best_settings$mzwid, minsamp=ud_minsamp, max=ud_max)
  
  message('about two')
  grouped = retcor.obiwarp(grouped, plottype="deviation", profStep=optimized_retcor_p$best_settings$profStep, gapInit=optimized_retcor_p$best_settings$gapInit, gapExtend=optimized_retcor_p$best_settings$gapExtend)
  
  message('about three')
  grouped = group(grouped,bw=optimized_retcor_p$best_settings$bw, minfrac=ud_minfrac, mzwid=optimized_retcor_p$best_settings$mzwid, minsamp=ud_minsamp, max=ud_max) # re-group data after RT alignment
  
  message('about to run fillPeaks()')
  results_chrom = fillPeaks(grouped, method='chrom') # ask which methods should be used
  #results_basic_batchI_MSW = fillPeaks(grouped, method='MSW') # this one gives an error.
  
  message('----')
}

setwd(results_dir)
save(xcms_optimized_parameters, file = "xcms_output.RData")
