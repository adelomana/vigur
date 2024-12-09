#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("ReactomePA")
#BiocManager::install("tictoc")

#
# 0. load libraries
#
library(crayon)
library(clusterProfiler)
library(enrichplot)
library(tictoc)
library(viridis)
library(ggplot2)

#
# 2. read files and generate lists of genes
#
filename = 'colored_proteins_04.tsv' 
df = read.csv(filename, sep='\t', header=TRUE)
subset = df[df$Color_in_plot == 'blue', ]
ensemblIDs = unique(unlist(strsplit(subset$ENSEMBL, ",")))
print(dim(df))
print(dim(subset))
print(length(ensemblIDs))
deps04blue = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
dim(deps04blue)

filename = 'colored_proteins_04.tsv' 
df = read.csv(filename, sep='\t', header=TRUE)
subset = df[df$Color_in_plot == 'red', ]
ensemblIDs = unique(unlist(strsplit(subset$ENSEMBL, ",")))
print(dim(df))
print(dim(subset))
print(length(ensemblIDs))
deps04red = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
dim(deps04red)

filename = 'colored_proteins_24.tsv' 
df = read.csv(filename, sep='\t', header=TRUE)
subset = df[df$Color_in_plot == 'blue', ]
ensemblIDs = unique(unlist(strsplit(subset$ENSEMBL, ",")))
print(dim(df))
print(dim(subset))
print(length(ensemblIDs))
deps24blue = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
dim(deps24blue)

filename = 'colored_proteins_24.tsv' 
df = read.csv(filename, sep='\t', header=TRUE)
subset = df[df$Color_in_plot == 'red', ]
ensemblIDs = unique(unlist(strsplit(subset$ENSEMBL, ",")))
print(dim(df))
print(dim(subset))
print(length(ensemblIDs))
deps24red = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
dim(deps24red)

geneLists = list('DEPs t = 4 h down'=deps04blue$ENTREZID, 
                 'DEPs t = 4 h up'=deps04red$ENTREZID, 
                 'DEPs t = 24 h down'=deps24blue$ENTREZID, 
                 'DEPs t = 24 h up'=deps24red$ENTREZID 
                )

geneLists = list(
                 'DEPs t = 24 h down'=deps24blue$ENTREZID, 
                 'DEPs t = 24 h up'=deps24red$ENTREZID 
)

#
# 3. run the analysis on different Ontologies
#
# this step takes surprisingly long time---
significance_threshold = 0.05
tic()
#ck = compareCluster(geneLists, fun="enrichPathway", pvalueCutoff=significance_threshold)
#ck = compareCluster(geneLists, fun='enrichGO', OrgDb='org.Hs.eg.db', pvalueCutoff=significance_threshold)
ck = compareCluster(geneLists, fun='enrichKEGG', pvalueCutoff=significance_threshold)
toc()
#cks = simplify(ck)
cks=ck

p1 = dotplot(cks, size='count', showCategory=10, font.size=8) 
print(p1)

# and I have a preference for cividis, but this is just personal preference
my_log_breaks = round(log10(significance_threshold)):round(log10(min(cks@compareClusterResult$p.adjust)))
my_log_breaks = my_log_breaks[ my_log_breaks %% 3 == 0] 
my_breaks = 10**my_log_breaks
p5 = p1 +  scale_fill_viridis(direction=-1, trans="log", breaks=my_breaks, option='cividis')
print(p5)
ggsave('dotplot.png', dpi=600)
dev.off()
# arguably this plot communicates best data patterns, IMHO

# importantly, store your fuctional enrichment in a form of table which will be a supplementary file of your paper
storage_file = 'clusterProfiler_enrichments.tsv'
write.table(cks@compareClusterResult, storage_file, quote=FALSE, sep='\t')
