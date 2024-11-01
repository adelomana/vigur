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
filename = 'deps.04.blue.tsv'
df = read.table(filename, sep='\t', header=TRUE)
ensemblIDs = df$ENSEMBL
deps04blue = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
dim(deps04blue)

filename = 'deps.04.red.tsv'
df = read.table(filename, sep='\t', header=TRUE)
ensemblIDs = df$ENSEMBL
deps04red = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
dim(deps04red)

filename = 'deps.24.blue.tsv'
df = read.table(filename, sep='\t', header=TRUE)
ensemblIDs = df$ENSEMBL
deps24blue = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
dim(deps24blue)

filename = 'deps.24.red.tsv'
df = read.table(filename, sep='\t', header=TRUE)
ensemblIDs = df$ENSEMBL
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
ck = compareCluster(geneLists, 
                    fun="enrichPathway", 
                    pvalueCutoff=significance_threshold)
toc()

p1 = dotplot(ck, size='count', showCategory=20, font.size=6) 
print(p1)

# and I have a preference for cividis, but this is just personal preference
my_log_breaks = round(log10(significance_threshold)):round(log10(min(ck@compareClusterResult$p.adjust)))
my_log_breaks = my_log_breaks[ my_log_breaks %% 2 == 0] 
my_breaks = 10**my_log_breaks
p5 = p1 +  scale_fill_viridis(direction=-1, trans="log", breaks=my_breaks, option='cividis')
print(p5)
# arguably this plot communicates best data patterns, IMHO

# importantly, store your fuctional enrichment in a form of table which will be a supplementary file of your paper
storage_file = 'clusterProfiler_enrichments.tsv'
write.table(ck@compareClusterResult, storage_file, quote=FALSE, sep='\t')
