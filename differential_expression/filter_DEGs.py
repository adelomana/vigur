###
### This script filters DEGs following the rules: (1) absolute log2 fold-change greater than one compared to initial conditions; (2) expression greater than five TPM; (3) relative standard error of the mean across replicates lower than one third; and (4) P < 0.05 and P-adjusted < 0.1 as called by DESeq2 [1] version 1.26.0.
###

import sys, numpy

###
### FUNC
###

def read_DEGs(DEGs, experiment, concentration):

    '''
    Returns a dictionary of DEGs.
    '''

    for time in time_tags:
        DEGs[experiment][concentration][time] = {}
        working_file = DESeq2_folder + experiment + '_' + concentration + '_time_' + time + '_vs_zero.csv'
        with open(working_file, 'r') as f:
            next(f)
            for line in f:
                v = line.split(',')
                gene_name = v[0].replace('"', '')
                log2FC = float(v[2])
                if v[-2] == 'NA':
                    pvalue = 1
                else:
                    pvalue = float(v[-2])
                if v[-1] == 'NA\n':
                    adjusted = 1
                else:
                    adjusted = float(v[-1])

                if (pvalue <= 0.05) and (adjusted <= 0.1):
                    DEGs[experiment][concentration][time][gene_name] = (log2FC, pvalue, adjusted)
                    
        # printing
        a = len(list(DEGs[experiment][concentration][time].keys()))
        print('{} {} {} DEGs: {}'.format(experiment, concentration, time, a))
    print()

    return DEGs

###
### MAIN
###

# 0. use-defined variables
DESeq2_folder = '/Volumes/sand/vigur/results/deseq2/'
experiment_tags = ['experiment_two', 'experiment_three']
concentration_tags = ['concentration_zero', 'concentration_half', 'concentration_five', 'concentration_fifty']
time_tags = ['four', 'twentyfour']
expression_file = '/Volumes/sand/vigur/results/expression/experiment_both_expression.txt'

# 1.1. define the DEGs across experimental design
print('define DEGs')
DEGs = {}
for experiment in experiment_tags:
    DEGs[experiment] = {}
    for concentration in concentration_tags:
        DEGs[experiment][concentration] = {}
        DEGs = read_DEGs(DEGs, experiment, concentration)

# 1.2. define expression
print('define expression')

with open(expression_file, 'r') as f:
    first_line = f.readline()
    metadata = first_line.split('\t')[1:]

print(len(metadata))

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
