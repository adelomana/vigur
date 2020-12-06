### IPO2

library(IPO)


data_dir = '/home/adrian/projects/vigur/data/metabolomics/pools/basic/'
peak_set = xcmsSet(datadir,  
                   ppm = 59.9, 
                   prefilter = c(3,100),
                   #                                      method= 'centWave',peakwidth= c(5.25, 21), snthresh = 10, 
                   #                                      mzCenterFun = "wMean", integrate = 1, mzdiff = -0.0355, noise = 0, fitgauss = FALSE,
                   #                                      polarity = "positive", verbose.columns = F) 

# # First need to generate a xcms set object, use the optimized parameters:
# inpdirPOS <- c("./PooledPos")
# optimizedXcmsSetObjectPOS <- xcmsSet(inpdirPOS,  ppm = 59.9, prefilter = c(3,100),
#                                      method= 'centWave',peakwidth= c(5.25, 21), snthresh = 10, 
#                                      mzCenterFun = "wMean", integrate = 1, mzdiff = -0.0355, noise = 0, fitgauss = FALSE,
#                                      polarity = "positive", verbose.columns = F) 
# 
# 
# time.RetGroup <- system.time({ # measuring time
#   resultRetcorGroup <-
#     optimizeRetGroup(xset = optimizedXcmsSetObjectPOS, 
#                      params = retcorGroupParameters, 
#                      nSlaves = 1, 
#                      subdir = NULL,
#                      plot = TRUE)
# })
# 
# # Too see the best settings:
# resultRetcorGroup$best_settings
# 
# # Best settings for acid POS:
# # gapInit: 0.4
# # gapExtend: 2.7
# # profStep: 1
# # bw: 4
# # mzwid: 0.035
# # minfrac: 0.7 - note this is completely up to the user....