import pandas, seaborn
import matplotlib, matplotlib.pyplot

###
### 0. user defined variables
###

data_file = '/home/adrian/projects/vigur/results/deseq2_filtered/strict_union_experiment_three.tsv'

### 1. read
df = pandas.read_csv(data_file, sep='\t')
df.set_index('gene_name', drop=True, inplace=True)
print(df)

#seaborn.heatmap(df, annot=True)

matplotlib.pyplot.pcolor(df, cmap='bwr', vmin=-3.5, vmax=3.5)
#matplotlib.pyplot.imshow(df, cmap='bwr', vmin=-3.5, vmax=3.5)

matplotlib.pyplot.savefig('figure.pdf')
