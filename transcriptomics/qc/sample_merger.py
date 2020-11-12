import sys

###
### This script reads kallisto folders and creates a single matrix file for quality control purposes.
###

# 0. user-defined variables
kallisto_dir = '/Volumes/sand/vigur/data/kallisto_shared_folders/'
results_file = '/Volumes/sand/vigur/results/expression/experiment_both_expression.txt'
metadata_file = '/Volumes/sand/vigur/data/metadata/vigur_metadata_experiment_both.tsv'

# 1. build metadata annotation
rosetta = {}

with open(metadata_file, 'r') as f:
    next(f)
    for line in f:
        v = line.split('\t')
        if v[0] != '':
            sample_id = v[0]
            label = 'experiment_{}_concentration_{}_time_{}_replicate_{}'.format(v[1], v[3], v[2], v[4].replace('\n', ''))
            rosetta[sample_id] = label
sample_names = list(rosetta.keys())

# 2. read expression
expression = {}
transcript_names = []

# read feature names
input_file = '/Volumes/sand/vigur/data/kallisto_shared_folders/RSS_HLMV_1/abundance.tsv'
with open(input_file, 'r') as f:
    next(f)
    for line in f:
        v = line.split('\t')
        transcript_names.append(v[0])
print('{} transcripts found'.format(len(transcript_names)))

# read tpm
for sample in sample_names:
    expression[rosetta[sample]] = {}
    input_file = kallisto_dir + sample + '/abundance.tsv'
    with open(input_file, 'r') as f:
        next(f)
        for line in f:
            v = line.split('\t')
            transcript_name = v[0]
            tpm = v[-1].replace('\n', '')
            expression[rosetta[sample]][transcript_name] = tpm

# 3. write file
with open(results_file, 'w') as f:

    f.write('transcriptID')
    for sample in sample_names:
        f.write('\t{}'.format(rosetta[sample]))
    f.write('\n')

    for transcript_name in transcript_names:
        f.write('{}'.format(transcript_name))
        for sample in sample_names:
            f.write('\t{}'.format(expression[rosetta[sample]][transcript_name]))
        f.write('\n')
        
