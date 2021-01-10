###
### This script check if any of the DEGs are dose responders
###

# read expression data and DEGs names from filter, so I know when, how and which trend they have
# compute their log10 round TPM + 1 against 0, 0.5, 5 and 50.
# define an R value as threshold

import pandas, sys, os, numpy
import scipy, scipy.stats

import pyensembl
annotation = pyensembl.EnsemblRelease(100)

import matplotlib, matplotlib.pyplot
matplotlib.rcParams.update({'font.size':20, 'font.family':'FreeSans', 'xtick.labelsize':20, 'ytick.labelsize':20})

###
### 0. user-defined variables
###

candidates_file = '/home/adrian/projects/vigur/data/transcriptomics/metabolic_candidates/Interesting redox and heparan sulphate proteins.csv'
tpm_file = '/home/adrian/projects/vigur/results/transcriptomics/deseq2/DESeq2_TPM_values.tsv'
metadata_file = '/home/adrian/projects/vigur/data/transcriptomics/metadata/vigur_metadata_experiment3.tsv'
DEGs_dir = '/home/adrian/projects/vigur/results/transcriptomics/deseq2_filtered/'
results_dir = '/home/adrian/projects/vigur/results/transcriptomics/dose/'

concentration_tags = ['zero', 'half', 'five', 'fifty']
time_tags = ['four', 'twentyfour']
trend_tags = ['up', 'down']
replicate_tags = ['A', 'B', 'C']

###
### 1. read data
###

### 1.1. create a datafrme of expression
expression = pandas.read_csv(tpm_file, sep='\t', index_col=0)

# convert expression dataframe columns into biological meaningful tags
metadata = pandas.read_csv(metadata_file, sep='\t')
for row in metadata.itertuples():
    old = row.sample
    new = 'time_{}_treatment_{}_replicate_{}'.format(row.time, row.treatment, row.replicate)
    new_names = {old:new}
    expression.rename(columns=new_names, inplace=True)

# read metabolic candidates
met_candidates = []
with open(candidates_file, 'r') as f:
    next(f)
    for line in f:
        v = line.split(',')
        met_candidates.append(v[1])

###
### 2. analysis
###

for time in time_tags:

    for trend in trend_tags:

        print(time, trend)

        ### 2.1. select all genes across concentrations
        working_genes = []
        for concentration in concentration_tags[1:]: # we don't want time zero for DEGs names
            reading_file = DEGs_dir + 'experiment_three_concentration_{}_time_{}_{}_filtered.tsv'.format(concentration, time, trend)
            if os.path.exists(reading_file):
                with open(reading_file, 'r') as f:
                    next(f)
                    next(f)
                    for line in f:
                        v = line.split('\t')
                        gene = v[0]
                        if gene not in working_genes:
                            working_genes.append(gene)
                print(concentration, len(working_genes))
        print()

        ### 2.2. iterate over the genes to see if there is enough correlation signal
        cases = []
        for concentration in concentration_tags:
            for replicate in replicate_tags:
                putative_case = 'time_{}_treatment_{}_replicate_{}'.format(time, concentration, replicate)
                if putative_case in expression.columns:
                    cases.append(putative_case)

        ### 2.3. slice dataframe by genes and cases
        reduced = expression[cases]
        print(reduced.shape)
        trimmed = reduced.loc[working_genes]
        print(trimmed.shape)


        ### 2.4. correlation analysis of expression
        association = pandas.DataFrame(index=working_genes, columns=['corr', 'delta'])
        for ensembl in working_genes:

            x = []; y = []
            # sort expression values depending on concentration
            for tag in trimmed.columns:

                tpm = trimmed.loc[ensembl][tag]
                rounded = numpy.around(tpm) + 1
                yvalue = numpy.log10(rounded)

                if 'treatment_zero' in tag:
                    xvalue = 0
                elif 'treatment_half' in tag:
                    xvalue = 1
                elif 'treatment_five' in tag:
                    xvalue = 2
                elif 'treatment_fifty' in tag:
                    xvalue = 3
                else:
                    raise ValueError("no concentration found")

                x.append(xvalue); y.append(yvalue)

            # linear regression analysis
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
            ypred = slope*numpy.array(x) + intercept
            delta = numpy.abs(numpy.max(ypred) - numpy.min(ypred))
            association.at[ensembl, 'corr'] = r_value
            association.at[ensembl, 'delta'] = delta

            # generate plot
            matplotlib.pyplot.plot(x, y, 'o', color='black', alpha=1/3, mew=0)
            matplotlib.pyplot.plot(x, ypred, '-', color='red', lw=2, alpha=2/3)
            title_tag = 'r = {:.2f}; P = {:.2e}'.format(r_value, p_value)
            matplotlib.pyplot.title(title_tag)
            matplotlib.pyplot.xticks([0, 1, 2, 3], ['0', '0.5', '5', '50'])
            matplotlib.pyplot.xlabel('Adrenaline')
            matplotlib.pyplot.ylabel('log10 TPM')
            figure_file = '{}time_{}_trend_{}_{}.pdf'.format(results_dir, time, trend, ensembl)
            matplotlib.pyplot.tight_layout()
            matplotlib.pyplot.savefig(figure_file)
            matplotlib.pyplot.close()

        ### 2.5. slice dataframe into groups
        strong = association[(numpy.abs(association['corr'].abs() > 0.8)) & (association['delta'] >= numpy.log10(2*2))]
        fair = association[(numpy.abs(association['corr'].abs() >= 0.8)) & (association['delta'] >= numpy.log10(2)) & (association['delta'] < numpy.log10(2*2))]
        weak = association[(numpy.abs(association['corr'].abs() < 0.8)) & (association['delta'] < numpy.log10(2*2))]

        ### 2.6. print dataframes for figure annotations
        gene_names = []
        for ensembl in strong.index:
            name = annotation.gene_name_of_gene_id(ensembl)
            gene_names.append(name)
            print(ensembl, name)
        print(gene_names)
        strong.insert(0, 'geneID', gene_names)

        gene_names = []
        for ensembl in fair.index:
            name = annotation.gene_name_of_gene_id(ensembl)
            gene_names.append(name)
        fair.insert(0, 'geneID', gene_names)

        print('strong')
        print(strong)
        print()
        print('fair')
        print(fair)
        print('weak')
        print(weak)
        print()

        ### 2.7. plot association results
        matplotlib.pyplot.plot(strong['corr'], strong['delta'], 'o', color='orange', alpha=2/3, mew=0, ms=8)
        matplotlib.pyplot.plot(fair['corr'], fair['delta'], 'o', color='green', alpha=1/2, mew=0, ms=6)
        matplotlib.pyplot.plot(weak['corr'], weak['delta'], 'o', color='black', alpha=1/3, mew=0, ms=4)

        for i in range(len(strong.index)):
            matplotlib.pyplot.text(strong['corr'][i], strong['delta'][i], strong['geneID'][i], fontsize=4)
        for i in range(len(fair.index)):
            matplotlib.pyplot.text(fair['corr'][i], fair['delta'][i], fair['geneID'][i], fontsize=4)

        if trend == 'up':
            matplotlib.pyplot.xlim([0.6, 1])
        if trend == 'down':
            matplotlib.pyplot.xlim([-1, -0.6])
        matplotlib.pyplot.ylim([0, 1.7])
        matplotlib.pyplot.xlabel('r value')
        matplotlib.pyplot.ylabel('delta')
        matplotlib.pyplot.grid(alpha=0.5)
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.savefig('{}association{}{}.svg'.format(results_dir, time, trend))
        matplotlib.pyplot.close()

        # overlap with metabolic model selected candidates (provided by Sarah)
        print('checking overlap...')
        print('candidates: {}; strong: {}'.format(len(met_candidates), len(strong)))
        print('candidates: {}; weak: {}'.format(len(met_candidates), len(weak)))
        a = list(set(list(strong.index)) & set(met_candidates))
        b = list(set(list(fair.index)) & set(met_candidates))
        c = list(set(list(weak.index)) & set(met_candidates))
        print('intersect', len(a), len(b), len(c))
        print(a, b, c)
        print()
