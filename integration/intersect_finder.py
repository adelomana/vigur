###
### 0. user-defined variables
###

rosetta = '/home/adrian/projects/vigur/results/integration/gProfiler_hsapiens_2-12-2021_11-19-41 AM.csv'


# 1 read the metabolic models in ensembl form
ensembl_IDs = []
with open(rosetta, 'r') as f:
    for line in f:
        v = line.split(',')
        field = v[1]
        value = field.replace('"', '')
        ensembl_IDs.append(value)
ensembl_IDs = list(set(ensembl_IDs))
ensembl_IDs.remove('nan')
print('{} converted IDs obtained'.format(len(ensembl_IDs)))
