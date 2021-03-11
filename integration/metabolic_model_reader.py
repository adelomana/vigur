###
### this script intersects transcriptional dose responders with metabolic model genes
###

import sys

###
### 0. user-defined variables
###

metabolic_genes_file = '/home/adrian/gd15/HI/research/vigur/data/metabolic_model_genes/HPMVECmodelgenes.csv'
dose_responders_file = 'xx'

###
### 1. read data
###

### 1.1. read metabolic genes
metabolic_model_ids = []
with open(metabolic_genes_file, 'r') as f:
    for line in f:
        v = line.split()
        ID = v[0]
        brokenID = ID.split('.')
        if len(brokenID) != 1:
            short = brokenID[0]
            metabolic_model_ids.append(short)

metabolic_model_ids = list(set(metabolic_model_ids))
print('found {} metabolic gene ids out of the original list of 2,006'.format(len(metabolic_model_ids)))

### 1.2. read dose responders

###
### 2. analysis
###

### 2.1. convert metabolic genes into ensembl ids
#metabolic_model_ensembls = []
#for gi in metabolic_model_ids:
    #print(gi)
