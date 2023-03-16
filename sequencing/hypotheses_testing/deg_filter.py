###
### This script filters DEGs following the rules: (1) absolute log2 fold-change greater than one compared to initial conditions; (2) expression greater than 2 TPM; (3) relative standard error of the mean across replicates lower than 1/2; and (4) P < 0.05 and P-adjusted < 0.1 as called by DESeq2 [1] version 1.26.0.
###

import sys, numpy

import matplotlib, matplotlib.pyplot
matplotlib.rcParams.update({'font.size':20, 'font.family':'FreeSans', 'xtick.labelsize':20, 'ytick.labelsize':20})

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
                time = v[1]
                treatment = v[2]
                replicate = v[3].replace('\n', '')

                metadata[sample] = ('experiment_three', 'concentration_'+treatment, 'time_'+time, 'replicate_'+replicate)

    return metadata

def read_DEGs(DEGs, experiment, concentration, time, trend):

    '''
    Returns a dictionary of DEGs.
    '''

    working_file = DESeq2_folder + experiment + '_' + concentration + '_' + time + '_' + trend + '.tsv'

    with open(working_file, 'r') as f:
        next(f)
        for line in f:
            v = line.split('\t')

            ensembl = v[0]
            biotype = v[7]
            description = v[8]
            hgnc_symbol = v[9].replace('\n', '')
            basemean = float(v[1])
            logFC = float(v[2])
            pvalue = float(v[5])
            adjusted = float(v[6])

            info = (ensembl, hgnc_symbol, biotype, description, basemean, logFC, pvalue, adjusted)

            DEGs[experiment][concentration][time][trend].append(info)

    print('{} {} {} {} \t DEGs: {}'.format(experiment, concentration, time, trend, len(DEGs[experiment][concentration][time][trend])))

    return DEGs

###
### MAIN
###

#
# 0. use-defined variables
#

metadata_file = '/home/adrian/projects/vigur/data/sequencing/metadata/metadata.tsv'
expression_file = '/home/adrian/projects/vigur/results/sequencing/tpm/DESeq2_TPM_values.tsv'
filtered_folder = '/home/adrian/projects/vigur/results/sequencing/DEGs_filtered/'
DESeq2_folder = '/home/adrian/projects/vigur/results/sequencing/DEGs/'

experiment_tags = ['run_72', 'run_73']
concentration_tags = ['treatment_mix', 'treatment_five_epi', 'treatment_five_nor', 'treatment_ilo_only', 'treatment_mix_plus_ilo']
time_tags = ['time_four', 'time_twentyfour']
trend_tags = ['up', 'down']

expression_threshold = 2
discrete_fc_threshold = 1
noise_threshold = 1/2

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
                try:
                    DEGs = read_DEGs(DEGs, experiment, concentration, time, trend)
                except:
                    pass
sys.exit()

# 1.2. define metadata
print('define metadata')
metadata = metadata_reader()
sample_IDs = list(metadata.keys())

# 1.3. define expression
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
# 2. filter genes
#
print('filter options:')
print('\t expression threshold: {}'.format(expression_threshold))
print('\t discrete FC threshold: {}'.format(discrete_fc_threshold))
print('\t noise threshold: {}'.format(noise_threshold))

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
                    # define references depending on concentration labels
                    if concentration == 'concentration_zero':
                        if metadata[sample][0:3] == (experiment, 'concentration_zero', 'time_zero'):
                            reference_labels.append(sample)
                    else:
                        if metadata[sample][0:3] == (experiment, 'concentration_zero', time):
                            reference_labels.append(sample)

                print('INFO: experiment {}, concentration {}, time {}, trend {}'.format(experiment, concentration, time, trend))
                print('INFO: sample labels {}, reference labels {}'.format(sample_labels, reference_labels))

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

                    # filter 3: identify noisy genes
                    ref_int = numpy.around(ref) + 1
                    sam_int = numpy.around(sam) + 1
                    sem_ref = numpy.std(ref_int) / numpy.sqrt(len(ref_int))
                    rsem_ref = sem_ref / numpy.mean(ref_int)
                    sem_sam = numpy.std(sam_int) / numpy.sqrt(len(sam_int))
                    rsem_sam = sem_sam / numpy.mean(sam_int)
                    noise = numpy.max([rsem_ref, rsem_sam])

                    # filter 4: relevant expression change difference---at least 2 TPM or > 5% of change
                    expression_difference = numpy.abs(s - r)
                    five_per_cent_change = numpy.median([r, s]) * 0.05
                    relevance_threshold = numpy.max([five_per_cent_change, 2])

                    # selection
                    if abs_log2FC < discrete_fc_threshold:
                        including = False
                        info = '\t WARNING: small change gene discarded. Expression changes from {:.3f} ({}) to {:.3f} ({}), resulting in abs_log2FC {:.3f}. {}, {}'.format(r, den, s, num, abs_log2FC, case[1], case[3])
                        print(info)

                    if (including == True) and (top < expression_threshold):
                        including = False
                        info = '\t WARNING: low-expression gene discarded. Expression changes from {:.3f} to {:.3f}. {}, {}'.format(r, s, case[1], case[3])
                        print(info)

                    if (including == True) and (noise > noise_threshold):
                        including = False
                        info = '\t WARNING: noisy gene. Ref: {}, RSEM {:.3f}; Sam: {}, RSEM {:.3f}. {}, {}'.format(ref, rsem_ref, sam, rsem_sam, case[1], case[3])
                        print(info)

                    #if (including == True) and (expression_difference < relevance_threshold):
                    #    including = False
                    #    info = '\t WARNING: non-relevant expression difference for {}. Expression difference: {:.2f}. Reference {:.3f} and sample {:.3f}.'.format(case[1], expression_difference, r, s)
                    #    print(info)

                    if including == True:
                        content = list(case)
                        content.append(r); content.append(s); content.append(abs_log2FC)
                        container.append(content)

                # info about low-expression cases
                after = len(container)
                print('{} {} {} {} | DEGs final reduction from {} to {}\n'.format(experiment, concentration, time, trend, before, after))

                # store reduced version
                DEGs[experiment][concentration][time][trend] = container

#
# 3. write results
#

# 3.2. flag DEGs of zero concentration
for experiment in experiment_tags:
    union = []
    for time in time_tags:
        for trend in trend_tags:
            for concentration in concentration_tags:
                if concentration == 'concentration_zero':
                    time_markers = []
                    for i in range(len(DEGs[experiment][concentration][time][trend])):
                        time_markers.append(DEGs[experiment][concentration][time][trend][i][0])
                        DEGs[experiment][concentration][time][trend][i].append('yes')
                    print('time markers for {} {} {} {} | {}'.format(experiment, time, trend, concentration, len(time_markers)))
                else:
                    confounding_time_markers = 0
                    for i in range(len(DEGs[experiment][concentration][time][trend])):
                        if DEGs[experiment][concentration][time][trend][i][0] in time_markers:
                            DEGs[experiment][concentration][time][trend][i].append('yes')
                            confounding_time_markers = confounding_time_markers + 1
                        else:
                            DEGs[experiment][concentration][time][trend][i].append('no')
                    specific_DEGs = len(DEGs[experiment][concentration][time][trend]) - confounding_time_markers
                    print('\t {} | found confounding time markers: {}/{} out of {} DEGs, {} specific.'.format(concentration, confounding_time_markers, len(time_markers), len(DEGs[experiment][concentration][time][trend]), specific_DEGs))

# 3.2. store
for experiment in experiment_tags:
    union = []
    for concentration in concentration_tags:
        for time in time_tags:
            for trend in trend_tags:
                storage = filtered_folder + '{}_{}_{}_{}_filtered.tsv'.format(experiment, concentration, time, trend)
                with open(storage, 'w') as f:
                    f.write('{}_{}_{}_{}\n'.format(experiment, concentration, time, trend))
                    f.write('ENSEMBL\tGene name\tBiotype\tDescription\tBase mean\tlog2FC\tP value\tAdjusted P-value\tReference expression (TPM)\tSample expression (TPM)\tDiscrete abs(log2FC)\tConsistency across experiments in this particular comparison\tTime marker\n')
                    for content in DEGs[experiment][concentration][time][trend]:
                        line = ''
                        for element in content:
                            if isinstance(element, str) == False:
                                if numpy.abs(element) < 0.05 and element != 0.:
                                    sub = '{:.4E}'.format(element)
                                else:
                                    sub = '{:.4f}'.format(element)
                            else:
                                sub = element
                            line=line+sub+'\t'
                        line=line+'\n'
                        f.write(line)

                        if content[0] not in union:
                            union.append(content[0])

    # 3.2.1. store union of DEGs for heatmap using external tools, like R
    union_storage = filtered_folder + 'strict_union_log2TPMplusOne.tsv'
    print('DEG union for experiment {}: {}'.format(experiment, len(union)))
    experiment = 'experiment_three'

    # build a dictionary ENSEMBL - GENENAME for DEGs
    rosetta = {}
    for concentration in concentration_tags:
        for time in time_tags:
            for trend in trend_tags:
                for case in DEGs[experiment][concentration][time][trend]:
                    if case[0] in union:
                        rosetta[case[0]] = case[1]

    with open(union_storage, 'w') as g:

        # header
        g.write('ENSEMBL\tGene name')

        working_conditions = []
        working_conditions.append(['concentration_zero', 'time_zero'])

        g.write('\t[Cat.] = 0 uM | Time = 0 h')

        for concentration in concentration_tags:
            for time in time_tags:
                a = concentration.split('_')[-1]
                if a == 'zero':
                    c = '0'
                elif a == 'half':
                    c = '0.5'
                elif a == 'five':
                    c = '5'
                elif a == 'fifty':
                    c = '50'
                if time  == 'time_four':
                    b = '4'
                else:
                    b = '24'
                full_tag = '[Cat.] = {} uM | Time = {} h'.format(c, b)
                working_conditions.append([concentration, time])
                g.write('\t{}'.format(full_tag))
        g.write('\n')

        # each line
        for ensemblID in union:
            g.write('{}'.format(ensemblID))

            # add the gene name
            g.write('\t{}'.format(rosetta[ensemblID]))

            for working_condition in working_conditions:
                concentration = working_condition[0]
                time = working_condition[1]

                # define sample labels
                sample_labels = []
                for sample in sample_IDs:
                    if metadata[sample][0:3] == (experiment, concentration, time):
                        sample_labels.append(sample)

                # gather TPM expression
                sam = [expression[label][ensemblID] for label in sample_labels]

                # compute value
                s = numpy.median(sam)
                value = numpy.round(numpy.log2(s+1), 3)
                #value = int(numpy.around(s))

                # write value
                g.write('\t{}'.format(value))

            # add a new line
            g.write('\n')
    g.close()

    # 3.2.1. store union of DEGs as log2FC for external tools
    union_storage = filtered_folder + 'strict_union_log2FC.tsv'
    with open(union_storage, 'w') as g:

        # header
        g.write('ENSEMBL\tGene name')

        for concentration in concentration_tags[1:]:
            for time in time_tags:

                if concentration == 'concentration_half':
                    formata = '0.5'
                elif concentration  == 'concentration_five':
                    formata = '5'
                elif concentration == 'concentration_fifty':
                    formata = '50'

                if time == 'time_four':
                    formatb = '4'
                elif time == 'time_twentyfour':
                    formatb = '24'
                word = '[Cat.] = {} | Time = {} h'.format(formata, formatb)
                g.write('\t{}'.format(word))
        g.write('\n')

        for ensemblID in union:

            # names
            g.write('{}'.format(ensemblID))
            g.write('\t{}'.format(rosetta[ensemblID]))

            for concentration in concentration_tags[1:]:
                for time in time_tags:

                    # define sample labels
                    sample_labels = []; reference_labels = []
                    for sample in sample_IDs:
                        if metadata[sample][0:3] == (experiment, concentration, time):
                            sample_labels.append(sample)
                        if metadata[sample][0:3] == (experiment, 'concentration_zero', time):
                            reference_labels.append(sample)

                    # gather TPM expression
                    ref = [expression[label][ensemblID] for label in reference_labels]
                    sam = [expression[label][ensemblID] for label in sample_labels]
                    r = numpy.median(ref); s = numpy.median(sam)

                    # log2 FC
                    num = numpy.around(s) + 1
                    den = numpy.around(r) + 1
                    fc = num/den
                    log2FC = numpy.round(numpy.log2(fc), 2)

                    # write value
                    g.write('\t{}'.format(log2FC))

            # add a new line
            g.write('\n')
