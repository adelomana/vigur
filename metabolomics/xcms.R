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

basic_batchI_dir = paste(data_dir, 'basic_batchI_dir', sep='')
basic_batchII_dir = paste(data_dir, 'basic_batchII_dir', sep='')
negative_batchI_dir = paste(data_dir, 'negative_batchI_dir', sep='')
negative_batchII_dir = paste(data_dir, 'negative_batchII_dir', sep='')
positive_batchI_dir = paste(data_dir, 'positive_batchI_dir', sep='')
positive_batchII_dir = paste(data_dir, 'positive_batchII_dir', sep='')

# 0.2. parameters
nThreads = 20

# basic_peak_width = c(3.0116056878683, 19.75823740928)
# basic_ppm = 205.562
# basic_mzdiff = -0.2342

basic_peak_width = c(8.074653184, 21.352)
basic_ppm = 84.8
basic_mzdiff = -0.06755

negative_peak_width = c(3, 20)
negative_ppm = 200
negative_mzdiff = -0.2

positive_peak_width = c(3, 20)
positive_ppm = 200
positive_mzdiff = -0.2

###
### 1. analysis
###

### 1.1. basic 
tic()
hlutur = xcmsSet(basic_batchI_dir,  
                 peakwidth=basic_peak_width,
                 ppm=basic_ppm,
                 mzdiff=basic_mzdiff, 
                 polarity='negative', # basic mode runs at negative polarity
                 BPPARAM=MulticoreParam(workers=nThreads),
                 #
                 prefilter=c(3, 100),
                 method='centWave',
                 snthresh=10, 
                 mzCenterFun="wMean", 
                 integrate=1,
                 noise=0, 
                 fitgauss=FALSE)

ud_minfrac = (6-1)/36 # at least 6 samples are biologically similar. 30 samples and 6 pools, 36 in total
ud_minfrac
ud_minsamp = 1
ud_max = 50

retcor_boundaries = getDefaultRetGroupStartingParams()
print('original boundaries')
print(retcor_boundaries)
retcor_boundaries$bw = c(0.1, 10) #  22-38 for HPLC, lower bounds for UPLC, which we are using. Freyrs bounds were (4, 25)
retcor_boundaries$minfrac = ud_minfrac
retcor_boundaries$mzwid = c(0.005, 0.2)
retcor_boundaries$profStep = c(0, 1)
retcor_boundaries$gapInit = c(0, 2)
#retcor_boundaries$gapExtend = c(1, 50)

optimized_retcor_p = optimizeRetGroup(xset=hlutur, params=retcor_boundaries, nSlaves=nThreads)

print('Optimization results:')
print(c('bw', retcor_boundaries$bw, optimized_retcor_p$best_settings$bw)) # (2,10) 2   
print(c('minfrac', retcor_boundaries$minfrac, optimized_retcor_p$best_settings$minfrac)) # (0.05, 0.5) 0.05
print(c('mzwid', retcor_boundaries$mzwid, optimized_retcor_p$best_settings$mzwid)) # 0.015 0.035 0.015
print(c('profStep', retcor_boundaries$profStep, optimized_retcor_p$best_settings$profStep)) # 0.5 1 0.5
print(c('gapInit', retcor_boundaries$gapInit, optimized_retcor_p$best_settings$gapInit)) # 0 0.8 0.8
print(c('gapExtend', retcor_boundaries$gapExtend, optimized_retcor_p$best_settings$gapExtend)) # 2.1 2.7 2.4

print(optimized_retcor_p$best_settings)

grouped = group(hlutur, bw=optimized_retcor_p$best_settings$bw, minfrac=ud_minfrac, mzwid=optimized_retcor_p$best_settings$mzwid, minsamp=ud_minsamp, max=ud_max)
grouped = retcor.obiwarp(grouped, plottype="deviation", profStep=optimized_retcor_p$best_settings$profStep, gapInit=optimized_retcor_p$best_settings$gapInit, gapExtend=optimized_retcor_p$best_settings$gapExtend)
grouped = group(grouped,bw=optimized_retcor_p$best_settings$bw, minfrac=ud_minfrac, mzwid=optimized_retcor_p$best_settings$mzwid, minsamp=ud_minsamp, max=ud_max) # re-group data after RT alignment
results_basic_batchI_chrom = fillPeaks(grouped, method='chrom') # ask which methods should be used
#results_basic_batchI_MSW = fillPeaks(grouped, method='MSW') # this one gives an error.
toc()
