library(biomaRt)
library(sleuth)
library(crayon) # so the messages are blue
library(tictoc)

#
# 0. user-defined variables
#
setwd("~/scratch/")

kallisto_dir = "/home/adrian/projects/vigur/results/sequencing/kallisto/kallisto.100"
metadata_file = '/home/adrian/projects/vigur/data/sequencing/metadata/metadata.tsv'
results_dir = '/home/adrian/projects/vigur/results/sequencing/DEGs_sleuth'

#
# 1. generate gene to transcript mapping
#
working_atributes = c('ensembl_transcript_id', 'ensembl_gene_id', 'gene_biotype', 'description', 'hgnc_symbol')
ensembl96 = useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", version=96)
table96 = getBM(attributes=working_atributes, mart=ensembl96)
colnames(table96)[1] = 'target_id'
#! listAttributes(mart=ensembl96)
View(table96)




#
# 2. read metadata
#
metadata = read.table(metadata_file, header=TRUE, sep="\t")
rownames(metadata) = metadata$sampleID # this line is necessary for later local selection of metadata during hypothesis testing
colnames(metadata)[2] = 'sample'
View(metadata)

#
# 3. work on different contrasts
#
cat(blue('quantify the effect of siMITF, IFN and their interaction'), fill=TRUE)

# 3.1. prepare metadata including paths
rules = metadata$run == '72' & metadata$experiment == 1 & metadata$treatment == 'zero' & metadata$time != 'twentyfour'
s2c = metadata[rules, ]
paths = file.path(kallisto_dir, s2c$sample, "abundance.h5")
s2c$path = paths
View(s2c)

# 3.2. prepare object for sleuth and plotting
cat(blue('preparing sleuth object...'), fill=TRUE)

tic()
so = sleuth_prep(s2c, 
                 target_mapping=table96, 
                 aggregation_column='ensembl_gene_id', 
                 read_bootstrap_tpm=TRUE, 
                 extra_bootstrap_summary=TRUE, 
                 read_bootstrap_tpm=TRUE, 
                 gene_mode=TRUE)
toc()

# 3.3. fit all models 
cat(blue('building models...'), fill=TRUE)
tic()
so = sleuth_fit(so, ~time, 'time')
so = sleuth_fit(so, ~1, 'reduced')
toc()