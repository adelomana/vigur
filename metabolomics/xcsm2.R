library('xcms')

# automatically detect the min frac


data_dir = '/home/adrian/projects/vigur/data/metabolomics/samples/pos_batch_I'

# variables
nThreads = 2

hlutur = xcmsSet(data_dir,  
                 peakwidth= c(4.536, 14.27),
                 ppm = 47.75,
                 mzdiff = -0.02465, 
                 polarity = "positive",
                 BPPARAM=MulticoreParam(workers=nThreads),
                 #
                 prefilter = c(3,100),
                 method= 'centWave',
                 snthresh = 10, 
                 mzCenterFun = "wMean", 
                 integrate = 1,
                 noise = 0, 
                 fitgauss = FALSE)