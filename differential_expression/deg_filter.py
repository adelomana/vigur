###
### This script filters DEGs following the rules: (1) absolute log2 fold-change greater than one compared to initial conditions; (2) expression greater than 2 TPM; (3) relative standard error of the mean across replicates lower than 1/2; and (4) P < 0.05 and P-adjusted < 0.1 as called by DESeq2 [1] version 1.26.0.
###

import sys, numpy, matplotlib_venn

import matplotlib, matplotlib.pyplot
matplotlib.rcParams.update({'font.size':20, 'font.family':'FreeSans', 'xtick.labelsize':20, 'ytick.labelsize':20})
#matplotlib.rcParams['pdf.fonttype']=42

###
### FUNCS
###

def extreme_responders_analysis(time, concentrations, figure_file):

    # f.1. iterate over concentrations for each time point compute x and y
    allx = []
    ally = []
    all_names = []

    for trend in trend_tags:
        enrichment_names = []; enrichment_IDs = []
        for concentration in concentrations:
            print(experiment, concentration, time, trend)

            for content in DEGs[experiment][concentration][time][trend]:
                # discrete log2FC considering sign
                if content[5] > 0:
                    x = content[10]
                elif content[5] < 0:
                    x = -content[10]
                else:
                    raise('DEG should not have a log2FC of zero')
                # log10 average expression
                y = numpy.log10(numpy.mean([content[8], content[9]]))
                allx.append(x)
                ally.append(y)
                all_names.append(content[1])
                enrichment_IDs.append(content[0])
                enrichment_names.append(content[1])
        enrichment_IDs = list(set(enrichment_IDs))
        enrichment_names = list(set(enrichment_names))
        print('\t Set of genes found: {}'.format(", ".join(enrichment_names)))
        for element in enrichment_IDs:
            print(element)

    ## f.2. calculate ranks of extreme outliers
    extreme_threshold = int(extreme_threshold_ratio * len(allx))
    print('searching for {} extreme responders which corresponds to a ratio of {} out of {} signficant events'.format(extreme_threshold, extreme_threshold_ratio, len(allx)))

    normalizedx = allx / numpy.max(numpy.abs(allx))
    normalizedy = ally / numpy.max(ally)
    objective = [numpy.sqrt(normalizedx[i]**2 + normalizedy[i]**2) for i in range(len(normalizedx))]
    positions = numpy.flip(numpy.argsort(objective))
    extreme_positions = positions[:extreme_threshold]

    epxo = []; epyo = []; epxg = []; epyg = []
    epo_names = []; epg_names =[]
    for pos in extreme_positions:
        if allx[pos] > 0:
             epxo = [allx[pos] for pos in extreme_positions if allx[pos] > 0]
             epyo = [ally[pos] for pos in extreme_positions if allx[pos] > 0]
             epo_names = [all_names[pos] for pos in extreme_positions if allx[pos] > 0]
        else:
            epxg = [allx[pos] for pos in extreme_positions if allx[pos] < 0]
            epyg = [ally[pos] for pos in extreme_positions if allx[pos] < 0]
            epg_names = [all_names[pos] for pos in extreme_positions if allx[pos] < 0]
        print(all_names[pos], allx[pos], ally[pos])
    greyx = []; greyy = []
    for i in range(len(allx)):
        if i not in extreme_positions:
            greyx.append(allx[i]); greyy.append(ally[i])

    ## f.3. plot
    matplotlib.pyplot.figure(figsize=(8, 8))
    matplotlib.pyplot.plot(greyx, greyy, 'o', color='black', alpha=1/8, mew=0, ms=10)
    matplotlib.pyplot.plot(epxo, epyo, 'o', color='orange', alpha=1/2, mew=0, ms=20)
    matplotlib.pyplot.plot(epxg, epyg, 'o', color='green', alpha=1/2, mew=0, ms=20)

    for i in range(len(epo_names)):
        matplotlib.pyplot.text(epxo[i], epyo[i], epo_names[i], horizontalalignment='left', verticalalignment='center')
    for i in range(len(epg_names)):
        matplotlib.pyplot.text(epxg[i], epyg[i], epg_names[i], horizontalalignment='left', verticalalignment='center')

    matplotlib.pyplot.xlabel('discrete log2 FC')
    matplotlib.pyplot.ylabel('log10 average expression')
    matplotlib.pyplot.tight_layout()
    #matplotlib.pyplot.savefig(figure_file)
    matplotlib.pyplot.savefig(figure_file, dpi=300)
    matplotlib.pyplot.clf()

    return None

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

project_dir = '/home/adrian/projects/vigur/'
DESeq2_folder = project_dir + 'results/deseq2/'
metadata_file = project_dir +'data/metadata/vigur_metadata_experiment_both.tsv'
expression_file = DESeq2_folder + 'DESeq2_TPM_values.tsv'
filtered_folder = project_dir + 'results/deseq2_filtered/'
venn_folder = project_dir + 'results/deseq2_venn/'

experiment_tags = ['experiment_two', 'experiment_three']
concentration_tags = ['concentration_zero', 'concentration_half', 'concentration_five', 'concentration_fifty']
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

                    # filter 3: identify noisy genes
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

# 3.1. identify consistent results
for concentration in concentration_tags:
    for time in time_tags:
        for trend in trend_tags:
            for experiment in experiment_tags:
                current = experiment
                other = [element for element in experiment_tags if element != current][0]

                # add consistency
                current_cases = list(set([element[0] for element in DEGs[current][concentration][time][trend]]))
                other_cases = list(set([element[0] for element in DEGs[other][concentration][time][trend]]))

                for i in range(len(DEGs[current][concentration][time][trend])):
                    case = DEGs[current][concentration][time][trend][i]

                    if case[0] in other_cases:
                        DEGs[current][concentration][time][trend][i].append('yes')
                    else:
                        DEGs[current][concentration][time][trend][i].append('no')

                # make Venn diagram
                if current == 'experiment_two':
                    print('make Venn diagram')
                    frame = numpy.sqrt(len(current_cases) + len(other_cases))/5
                    if trend == 'up':
                        the_colors = ('orange', 'red')
                    else:
                        the_colors = ('green', 'blue')
                    matplotlib.pyplot.figure(figsize=(frame, frame))
                    matplotlib_venn.venn2([set(current_cases), set(other_cases)], set_labels = ('Experiment 2', 'Experiment 3'), set_colors=the_colors)
                    title_tag = '{}_{}_trend_{}'.format(concentration, time, trend)
                    matplotlib.pyplot.title(title_tag)
                    figure_file = '{}.svg'.format(venn_folder+title_tag)
                    matplotlib.pyplot.savefig(figure_file)
                    matplotlib.pyplot.clf()

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

    # store union of DEGs
    union_storage = filtered_folder + 'union_{}.tsv'.format(experiment)
    print('DEG union for experiment {}: {}'.format(experiment, len(union)))
    with open(union_storage, 'w') as g:
        for element in union:
            g.write('{}\n'.format(element))

# 4. make a figure with extreme responders
# extreme responders should be the 5% of log2 FC plus average expression.
# highlight the top lower ranks, summing both ranks of average expression and log2 FC.

extreme_threshold_ratio = 0.05
experiment = 'experiment_three'

#########################################

time = 'time_four'
concentrations = ['concentration_zero']
figure_file = filtered_folder + 'figures/extreme_responders_time_four_without_adrenaline.png'
extreme_responders_analysis(time, concentrations, figure_file)

time = 'time_twentyfour'
concentrations = ['concentration_zero']
figure_file = filtered_folder + 'figures/extreme_responders_time_twentyfour_without_adrenaline.png'
extreme_responders_analysis(time, concentrations, figure_file)

# ADRENALINE HALF ##########################################
time = 'time_four'
concentrations = ['concentration_half']
figure_file = filtered_folder + 'figures/extreme_responders_time_four_concentration_half.png'
extreme_responders_analysis(time, concentrations, figure_file)

time = 'time_twentyfour'
concentrations = ['concentration_half']
figure_file = filtered_folder + 'figures/extreme_responders_time_twentyfour_concentration_half.png'
extreme_responders_analysis(time, concentrations, figure_file)

# ADRENALINE FIVE ##########################################
time = 'time_four'
concentrations = ['concentration_five']
figure_file = filtered_folder + 'figures/extreme_responders_time_four_concentration_five.png'
extreme_responders_analysis(time, concentrations, figure_file)

time = 'time_twentyfour'
concentrations = ['concentration_five']
figure_file = filtered_folder + 'figures/extreme_responders_time_twentyfour_concentration_five.png'
extreme_responders_analysis(time, concentrations, figure_file)

# ADRENALINE FIFTY ##########################################
time = 'time_four'
concentrations = ['concentration_fifty']
figure_file = filtered_folder + 'figures/extreme_responders_time_four_concentration_fifty.png'
extreme_responders_analysis(time, concentrations, figure_file)

time = 'time_twentyfour'
concentrations = ['concentration_fifty']
figure_file = filtered_folder + 'figures/extreme_responders_time_twentyfour_concentration_fifty.png'
extreme_responders_analysis(time, concentrations, figure_file)
