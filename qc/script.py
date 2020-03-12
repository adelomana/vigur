import sys

# 0. user-defined variables
kallisto_dir = '/Volumes/sand/vigur/data/kallisto_shared_folders/'
results_dir = '/Volumes/sand/vigur/results/expression/'

# 1. build annotationi
rosetta = {}
index = 0
treatments = ['0.0', '0.5', '5', '50']
times = ['0', '4', '24']

# add time point zero
for k in range(3):
    index = index + 1
    folder_tag = 'RSS_HLMV_{}'.format(index)
    label =  'RSS_HLMV{}_time0_treatment0.0_replicate{}'.format(index,index)
    print(folder_tag,label)
    rosetta[folder_tag] = label

for i in range(2):
    for j in range(4):
        for k in range(3):
            index = index + 1
            folder_tag = 'RSS_HLMV_{}'.format(index)
            label =  'RSS_HLMV{}_time{}_treatment{}_replicate{}'.format(index, times[i+1], treatments[j], k+1)
            print(folder_tag,label)
            rosetta[folder_tag] = label

print(rosetta)

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

# read tpm
for sample in rosetta:
    expression[rosetta[sample]] = {}
    input_file = kallisto_dir + sample + '/abundance.tsv'
    with open(input_file, 'r') as f:
        next(f)
        for line in f:
            v = line.split('\t')
            transcript_name = v[0]
            tpm = v[-1].replace('\n', '')
            expression[rosetta[sample]][transcript_name] = tpm

# 3.
