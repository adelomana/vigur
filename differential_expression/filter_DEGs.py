###
### This script filters DEGs following the rules: (1) absolute log2 fold-change greater than one compared to initial conditions; (2) expression greater than five TPM; (3) relative standard error of the mean across replicates lower than one third; and (4) P < 0.05 and P-adjusted < 0.1 as called by DESeq2 [1] version 1.26.0.
###

import sys, numpy

###
### FUNC
###

def read_DEGs(DEGs, experiment, concentration, time, trend):

    '''
    Returns a dictionary of DEGs.
    '''

    working_file = DESeq2_folder + experiment + '_' + concentration + '_' + time + '_vs_zero' + '_' + trend + '.tsv'
    
    with open(working_file, 'r') as f:
        next(f)
        for line in f:
            v = line.split('\t')
            
            ensembl = v[0]
            gene_name = v[1]
            biotype = v[2]
            description = v[3]
            basemean = float(v[4])
            logFC = float(v[5])
            pvalue = float(v[6])
            adjusted = float(v[7])

            info = (ensembl, gene_name, biotype, description, basemean, logFC, pvalue, adjusted)

            DEGs[experiment][concentration][time][trend].append(info)
                    
    print('{} {} {} {} \t DEGs: {}'.format(experiment, concentration, time, trend, len(DEGs[experiment][concentration][time][trend])))

    return DEGs

###
### MAIN
###

# 0. use-defined variables
DESeq2_folder = '/Users/alomana/projects/vigur/results/deseq2/'
experiment_tags = ['experiment_two', 'experiment_three']
concentration_tags = ['concentration_zero', 'concentration_half', 'concentration_five', 'concentration_fifty']
time_tags = ['time_four', 'time_twentyfour']
trend_tags = ['up', 'down']
expression_file = DESeq2_folder + 'DESeq2_TPM_values.tsv'

# 1.1. define the DEGs across experimental design
print('define DEGs')
DEGs = {}

for experiment in experiment_tags:
    DEGs[experiment] = {}
    for concentration in concentration_tags:
        DEGs[experiment][concentration] = {}
        for time in time_tags:
            DEGs[experiment][concentration][time] = {}
            for trend in trend_tags:
                DEGs[experiment][concentration][time][trend] = []

                DEGs = read_DEGs(DEGs, experiment, concentration, time, trend)
                
# 1.2. define expression
print('define expression')

with open(expression_file, 'r') as f:
    first_line = f.readline()
    metadata = first_line.split('\t')[1:]
    print(metadata)
    print(len(metadata))

sys.exit()

expression = {}
with open(expression_file, 'r') as f:
    next(f)
    for line in f:
        v = line.split('\t')
        gene_name = v[0]
        v = [float(element) for element in v[1:]]

        print(gene_name)
        print(v)
        print(len(v))
        
        #for experiment in experiment_tags:
    #for concentration in concentration_tags:
    #    for time in time_tags:
        
      #  sys.exit()
    
sys.exit()

# 2.1. define DEGs after filtering
print('filter concentration_zero DEGs')
for experiment in experiment_tags:
    
    # 2.1. define concentration_zero DEGs
    zero_genes = []
    for time in time_tags:
        sub = list(DEGs[experiment]['concentration_zero'][time].keys())
        print(experiment, time, len(sub))
        for gene in sub:
            if gene not in zero_genes:
                zero_genes.append(gene)
        print(len(zero_genes))

    # 2.2. filter concentration_zero DEGs in other concentrations
    for concentration in concentration_tags[1:]:
        for time in time_tags:
            a = len(list(DEGs[experiment][concentration][time]))
            print('\t {} {} {} \t before {}'.format(experiment, concentration, time, a))
            
            for zero_gene in zero_genes:
                if zero_gene in DEGs[experiment][concentration][time].keys():
                    del DEGs[experiment][concentration][time][zero_gene]

            a = len(list(DEGs[experiment][concentration][time]))
            print('\t {} {} {} \t after {}'.format(experiment, concentration, time, a))
            print()
            

# 2.3. filter for abs log2FC > 1
print('filter log2 FC')
for experiment in experiment_tags:
    for concentration in concentration_tags:
        for time in time_tags:
            small_change_genes = []

            a = len(list(DEGs[experiment][concentration][time]))
            print('\t {} {} {} \t before {}'.format(experiment, concentration, time, a))
            
            for gene in list(DEGs[experiment][concentration][time].keys()):
                fc = numpy.abs(DEGs[experiment][concentration][time][gene][0])
                if fc < 1:
                    small_change_genes.append(gene)

            print('\t {} small change genes detected'.format(len(small_change_genes)))

            for small_change_gene in small_change_genes:
                del DEGs[experiment][concentration][time][small_change_gene]
                
            a = len(list(DEGs[experiment][concentration][time]))
            print('\t {} {} {} \t after {}'.format(experiment, concentration, time, a))
            print()

# 2.4. filter for at least 5 TPMs


# 2.5. filter out noisy genes: relative standard error of the mean across replicates lower than one third
