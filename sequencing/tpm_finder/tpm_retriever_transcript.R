#
# -1. packages installation
#
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# BiocManager::install("tictoc")

library(biomaRt)
library(sleuth)
library(crayon) # so the messages are blue
library(tictoc)
  
#
# 0. user-defined variables
#
setwd("~/scratch/")

kallisto_dir = "/home/adrian/projects/vigur/results/sequencing/kallisto/kallisto.100"
results_dir = '/home/adrian/projects/vigur/results/sequencing/tpm'

#
# 1. generate gene to transcript mapping
#

# annotation defined from sleuth walk through, https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html
mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host='https://uswest.ensembl.org')
t2g = biomaRt::getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart=mart)
t2g = dplyr::rename(t2g, target_id=ensembl_transcript_id, ens_gene=ensembl_gene_id, ext_gene=external_gene_name)

#
# 2. read metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
dirnames

files = file.path(dirnames, 'abundance.h5')
print(files)
sample = sapply(strsplit(files, split='/',fixed=TRUE), function(x) (x[10]))
metadata = data.frame(sample)
metadata

#
# 3. create a sleuth object
#
s2c = metadata
paths = file.path(kallisto_dir, s2c$sample, "abundance.h5")
s2c$path = paths
print(s2c)

# prepare object for sleuth
cat(blue('preparing sleuth object...'), fill=TRUE)
tic()
so = sleuth_prep(s2c, 
                 target_mapping=t2g, 
                 aggregation_column='ens_gene', 
                 read_bootstrap_tpm=TRUE, 
                 extra_bootstrap_summary=TRUE)
toc()

#
# 4. store TPMs
#
cat(blue('storing'), fill=TRUE)
tpm_table = sleuth_to_matrix(so, 'obs_norm', 'tpm')
View(tpm_table)

write.csv(tpm_table, file.path(results_dir, 'sleuth_TPM_transcript.csv'))
