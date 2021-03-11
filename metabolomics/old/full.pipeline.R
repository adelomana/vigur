#!/usr/bin/env Rscript

source('pipeline.library.R')

###
### 0. user-defined variables
###

setwd('/home/adrian/hub/vigur/metabolomics/')

# paths
pool_data_dir = '/home/adrian/projects/vigur/data/metabolomics/pools/'
pool_results_dir = '/home/adrian/projects/vigur/results/metabolomics/pools/'

# variables
nThreads = 20

###
### 1. read info
###
cases = list.files(pool_data_dir)
cases = cases[1:3]

###
### 2. analysis
###
ipo_optimized_params_fullset = data.frame()

for (case_index in 1:length(cases)) {
  
  label = cases[case_index]

  the_message = paste(now(), '\t', 'working with', label)
  message(the_message)
  
  ###
  ### 2.1. IPO
  ###
  the_message = paste(now(), '\t\t', 'starting IPO...')
  message(the_message)

  subdir = paste(pool_data_dir, label, sep='')
  pool_files = list.files(subdir, recursive=TRUE, full.names=TRUE) 
  print(pool_files)  
  
  opt_ipo_params = ipo_caller(label, pool_files, pool_results_dir, nThreads)
  ipo_optimized_params_fullset = rbind(ipo_optimized_params_fullset, opt_ipo_params)
  
  the_message = paste(now(), '\t\t', 'IPO completed')
  message(the_message)
  
  ###
  ### 2.2. xcms
  ###
  the_message = paste(now(), '\t\t', 'starting xcms...')
  message(the_message)
  
  result = xcms_caller()
  
  the_message = paste(now(), '\t\t', 'xcms completed')
  message(the_message)
  
}