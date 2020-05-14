###
### This script filters DEGs following the rules: (1) absolute log2 fold-change greater than one compared to initial conditions; (2) expression greater than five TPM; (3) relative standard error of the mean across replicates lower than one third; and (4) P < 0.05 and P-adjusted < 0.1 as called by DESeq2 [1] version 1.26.0.
###

import sys, numpy

###
### FUNCS
###

def metadata_reader():

    metadata = {}
    with open(metadata_file, 'r') as f:
        next(f)
        for line in f:
            v = line.split('\t')
            if v[0] != '':

                sample = v[0]
                experiment = v[1]
                time = v[2]
                treatment = v[3]
                replicate = v[4].replace('\n', '')

                metadata[sample] = ('experiment_'+experiment, 'concentration_'+treatment, 'time_'+time, 'replicate_'+replicate)
    
    return metadata

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

#
# 0. use-defined variables
#

DESeq2_folder = '/Users/alomana/projects/vigur/results/deseq2/unfiltered/'
metadata_file = '/Users/alomana/projects/vigur/data/metadata/vigur_metadata_experiment_both.tsv'
expression_file = DESeq2_folder + 'DESeq2_TPM_values.tsv'
filtered_folder = '/Users/alomana/projects/vigur/results/deseq2/filtered/'

experiment_tags = ['experiment_two', 'experiment_three']
concentration_tags = ['concentration_zero', 'concentration_half', 'concentration_five', 'concentration_fifty']
time_tags = ['time_four', 'time_twentyfour']
trend_tags = ['up', 'down']

expression_threshold = 3
discrete_fc_threshold = 1
noise_threshold = 1/3

#
# 1. read data
#

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

# 1.2. define metadata
print('define metadata')
metadata = metadata_reader()
sample_IDs = list(metadata.keys())
                
# 1.2. define expression
print('define expression')

expression = {}
for sample in sample_IDs:
    expression[sample] = {}
    
with open(expression_file, 'r') as f:
    next(f)
    for line in f:
        v = line.split('\t')
        
        gene_name = v[0]
        v = [float(element) for element in v[1:]]

        for i in range(len(v)):
            value = v[i]
            sampleID = sample_IDs[i]
            expression[sampleID][gene_name] = value

#
# 2. analysis
#

# 2.1. filter
print('filter DEGs')

for experiment in DEGs:
    for concentration in DEGs[experiment]:
        for time in DEGs[experiment][concentration]:
            for trend in DEGs[experiment][concentration][time]:                

                # define sample and reference sample labels
                sample_labels = []; reference_labels = []
                for sample in sample_IDs:
                    if metadata[sample][0:3] == (experiment, concentration, time):
                        sample_labels.append(sample)
                    if metadata[sample][0:3] == (experiment, 'concentration_zero', 'time_zero'):
                        reference_labels.append(sample)

                # filters
                container = []
                before = len(DEGs[experiment][concentration][time][trend])
                for case in DEGs[experiment][concentration][time][trend]:
                    including = True
                    ensembl = case[0]

                    # gather TPM expression
                    ref = [expression[label][ensembl] for label in reference_labels]
                    sam = [expression[label][ensembl] for label in sample_labels]

                    # filter 1: identify low-expressed genes
                    r = numpy.median(ref); s = numpy.median(sam)
                    top = numpy.max([r, s])

                    # filter 2: identify fold-changes using discrete values
                    ###
                    ###            [round(x, epsilon)/epsilon ] + 1 
                    ###  FC = abs  -------------------------------- > 1
                    ###            [round(y, epsilon)/epsilon ] + 1
                    ###
                    ###
                    ###  epsilon = 1
                    num = numpy.around(s) + 1
                    den = numpy.around(r) + 1
                    fc = num/den
                    abs_log2FC = numpy.abs(numpy.log2(fc))

                    # filter 3: noisy genes, rsem > 1/3
                    ref_int = numpy.around(ref) + 1
                    sam_int = numpy.around(sam) + 1
                    sem_ref = numpy.std(ref_int) / numpy.sqrt(len(ref_int))
                    rsem_ref = sem_ref / numpy.mean(ref_int)
                    sem_sam = numpy.std(sam_int) / numpy.sqrt(len(sam_int))
                    rsem_sam = sem_sam / numpy.mean(sam_int)
                    noise = numpy.max([rsem_ref, rsem_sam])
                    
                    # selection
                    if abs_log2FC < discrete_fc_threshold:
                        including = False
                        info = 'WARNING: small change gene discarded.\nExpression changes from {:.3f} ({}) to {:.3f} ({}), resulting in abs_log2FC {:.3f}.\n{}, {}\n'.format(r, den, s, num, abs_log2FC, case[1], case[3])
                        print(info)
                        
                    if (including == True) and (top < expression_threshold):
                        including = False
                        info = 'WARNING: low-expression gene discarded.\nExpression changes from {:.3f} to {:.3f}.\n{}, {}\n'.format(r, s, case[1], case[3])
                        print(info)
                        
                    if (including == True) and (noise > noise_threshold):
                        including = False
                        info = 'WARNING: noisy gene.\nRef: {}, RSEM {:.3f}\nSam: {}, RSEM {:.3f} \n{}, {}\n'.format(ref, rsem_ref, sam, rsem_sam, case[1], case[3])
                        print(info)
                        
                    if including == True:
                        content = list(case)
                        content.append(r); content.append(s); content.append(abs_log2FC)
                        container.append(content)

                # info about low-expression cases
                after = len(container)
                print('{} {} {} {} | DEGs reduction from {} to {}'.format(experiment, concentration, time, trend, before, after))
                
                # write results
                storage = filtered_folder + '{}_{}_{}_{}_filtered.tsv'.format(experiment, concentration, time, trend)
                with open(storage, 'w') as f:
                    f.write('#ENSEMBL\tGene name\tBiotype\tDescription\tBase mean\tlog2FC\tP value\tAdjusted P-value\tReference expression (TPM)\tSample expression (TPM)\n')
                    for content in container:
                        info = '\t'.join([str(element) for element in content])
                        line = '{}\n'.format(info)
                        f.write(line)
