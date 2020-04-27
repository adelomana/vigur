
# 0. install sleuth
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("rhdf5", "devtools"))
devtools::install_github("pachterlab/sleuth")

# 1. load sleuth
library("sleuth")
