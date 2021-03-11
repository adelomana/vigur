ipo_caller = function(label, pool_files, pool_results_dir, nThreads){
  
  # just for testing
  pool_files = pool_files[1]
  
  originalwd = getwd()
  results_subdir = paste(pool_results_dir, label, sep='')
  dir.create(results_subdir)
  setwd(results_subdir)
  
  # explore parameters
  ppParameters = getDefaultXcmsSetStartingParams('centWave')
  ppParameters$min_peakwidth = c(2, 10)
  print(ppParameters)
  
  # freyrs boundaries
  #ppParameters$min_peakwidth = c(2, 7)
  #ppParameters$max_peakwidth = c(10, 20)
  #ppParameters$ppm = c(5, 50)
  #ppParameters$mzdiff = c(-0.05, 0.010)
  
  resultParameters = optimizeXcmsSet(files=pool_files, params=ppParameters, plot=TRUE, BPPARAM=MulticoreParam(workers=nThreads))
  
  # write results
  results_file = paste(label, '.tsv', sep='')
  df = as.data.frame(resultParameters$best_settings$parameters)
  write.table(df, file=results_file, sep='\t', quote=FALSE)
  
  setwd(originalwd)
  
  print(df)
  
  return(df)
}

xcms_caller = function(x){
  return('done')
}