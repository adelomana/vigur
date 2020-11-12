install.packages("diceR")
install.packages('pander')


library(diceR)
library(dplyr)
library(ggplot2)
library(pander)
data(hgsc)
hgsc <- hgsc[1:100, 1:50]

CC <- consensus_cluster(hgsc, nk = 3:4, p.item = 0.8, reps = 5,
                        algorithms = c("hc", "pam", "diana"))


co <- capture.output(str(CC))
strwrap(co, width = 80)


CC <- apply(CC, 2:4, impute_knn, data = hgsc, seed = 1)
CC_imputed <- impute_missing(CC, hgsc, nk = 4)
sum(is.na(CC))
#> [1] 5
sum(is.na(CC_imputed))


pam.4 <- CC[, , "PAM_Euclidean", "4", drop = FALSE]
cm <- consensus_matrix(pam.4)
dim(cm)
#> [1] 100 100

hm <- graph_heatmap(pam.4)
