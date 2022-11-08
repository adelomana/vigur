#
# 0. load libraries
#
library(crayon) # so messages can be in blue font
library(clusterProfiler)
library(enrichplot)

#
# 1. user-defined variables
#
gene_sets_dir = '/home/adrian/projects/vigur/results/proteomics'
output_dir = '/home/adrian/projects/vigur/results/proteomics/enrichment'

working_labels = list.files(gene_sets_dir, pattern='.txt')
working_labels

#
# 2. iterate gene sets
#
all_identifiers = list()
for (working_label in working_labels) {
  
  ## 2.1. retrieve genes
  working_file = paste(gene_sets_dir, working_label, sep='/')
  df = read.csv(working_file)
  identifiers = df$entrez_id
  cat(blue(paste('working with', working_label, 'with', length(identifiers), 'genes', sep=' ')), fill=TRUE)
  
  all_identifiers[[working_label]] = identifiers
  
  cat('\n')
}

#
# 3. run the analysis on different Ontologies
#

ck = compareCluster(all_identifiers, fun="enrichPathway", pvalueCutoff=0.05)
sm = pairwise_termsim(ck)
dotplot(ck, size='count', showCategory=10, title='enrichPathway')
emapplot(sm)
storage = paste(output_dir, 'clusterProfiler_enrichments_RP.csv', sep='/')
write.csv(ck@compareClusterResult, storage)

ck = compareCluster(all_identifiers, fun="enrichGO", pvalueCutoff=0.05, OrgDb='org.Hs.eg.db')
sm = pairwise_termsim(ck)
dotplot(ck, size='count', showCategory=10, title='enrichGO')
emapplot(sm)

ck = compareCluster(all_identifiers, fun="enrichKEGG", pvalueCutoff=0.05)
sm = pairwise_termsim(ck)
dotplot(ck, size='count', showCategory=10, title='enrichKEGG')
emapplot(sm)

ck = compareCluster(all_identifiers, fun="enrichDO", pvalueCutoff=0.05)
sm = pairwise_termsim(ck)
dotplot(ck, size='count', showCategory=5, title='enrichDO')
emapplot(sm)
