#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("heatmap3")
#BiocManager::install("gplots")

library(RColorBrewer)
library(pheatmap)

###
### 0. user-defined variables
###

data_file = '/home/adrian/projects/vigur/results/deseq2_filtered/strict_union_experiment_three.tsv'

###
### 1. read data
###

df = read.csv(data_file, sep='\t', row.names=1)
head(df)

###
### 2. analysis
###

dim(df)
df = subset(df, select = -c(zero4, zero24))
df = subset(df, select = -c(half24, five24, fifty24))

excluded = c()
for (gene in 1:nrow(df)) {
  if (max(abs(df[gene, ])) < 0.5){excluded = c(excluded, gene)}
}
df = df[-excluded, ]
dim (df)

data = as.matrix(df)

red_breaks = seq(0.5, 3.5, 0.1)
blue_breaks = seq(-3.5, -0.5, 0.1)
white_breaks = seq(-0.5, 0.5, 0.1)
white_breaks = white_breaks[-c(1, length(white_breaks))]
all_breaks = c(blue_breaks, white_breaks, red_breaks)
all_breaks

red_map = colorRampPalette(brewer.pal(9, "Reds"))(length(red_breaks))
blue_map = rev(colorRampPalette(brewer.pal(9, "Blues"))(length(red_breaks)))
white_map = rep('#FFFFFF', length(white_breaks)-1)
all_map = c(blue_map, white_map, red_map)
all_map

#pheatmap(data, col=all_map, breaks=all_breaks, clustering_distance_rows = "euclidean", clustering_method = "ward.D2")

###
### 4. consensus clustering
###

all_algorithms = c("diana", "pam")
cc = consensus_cluster(df, nk=4:8, p.item=0.9, reps=10, algorithms=all_algorithms, scale=FALSE)

pam.6 <- cc[, , "PAM_Euclidean", "6", drop = FALSE]
cm <- consensus_matrix(pam.6)
hm <- graph_heatmap(pam.6)
new_order = hm[[2]]$`PAM_Euclidean k=6`
new_order[, 'ensembl'] = rownames(df)
head(new_order, 20)
sorted = new_order[order(new_order$Cluster), ]
head(sorted, 20)

fin = df[match(sorted$ensembl, rownames(df)), ]
head(fin, 20)

pheatmap(fin, col=all_map, breaks=all_breaks, clustering_distance_rows = "euclidean", clustering_method = "ward.D2", cluster_cols=FALSE, cluster_rows=FALSE)

pheatmap(fin, col=all_map, breaks=all_breaks, clustering_distance_rows = "euclidean", clustering_method = "ward.D2")

