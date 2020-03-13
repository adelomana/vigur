import sys

###
### This script reads kallisto folders and creates a single matrix file for quality control purposes.
###

# 0. user-defined variables
kallisto_dir = '/Volumes/sand/vigur/data/kallisto_shared_folders/'
results_file = '/Volumes/sand/vigur/results/expression/experiment2.expression.txt'

# 1. build annotationi
rosetta = {}
index = 0
treatments = ['00.0', '00.5', '05.0', '50.0']
times = ['00', '04', '24']

# add time point zero
for k in range(3):
    index = index + 1
    folder_tag = 'RSS_HLMV_{}'.format(index)
    label =  'RSS_HLMV{}_time{}_treatment{}_replicate{}'.format(index, times[0], treatments[0], index)
    rosetta[folder_tag] = label

for i in range(2):
    for j in range(4):
        for k in range(3):
            index = index + 1
            folder_tag = 'RSS_HLMV_{}'.format(index)
            label =  'RSS_HLMV{}_time{}_treatment{}_replicate{}'.format(index, times[i+1], treatments[j], k+1)
            print(folder_tag,label)
            rosetta[folder_tag] = label

sample_names = list(rosetta.keys())
#sample_names.sort()

print('{} samples found'.format(len(sample_names)))

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
#transcript_names.sort()

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
        
