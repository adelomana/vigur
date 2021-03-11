###
### usage: time python kallisto_runner.py &> kallisto_runner_messages.txt
###

import sys, datetime, os

def kallisto_caller(working_labels):


    executable='time kallisto quant'
    options=' -i {} -o {} --bias -t {} --genomebam --gtf {} '.format(transcriptome_index, results_dir, threads, gtf_file)

    fastq_files = []
    for label in working_labels:
        for i in range(2):
            file = '{}{}1.R{}.fastq.gz'.format(clean_fastq_dir, label, i+1)
            fastq_files.append(file)

    fastq_files_string = ' '.join(fastq_files)

    command=executable+options+fastq_files_string

    print('')
    print(command)
    os.system(command)
    print('')

    return None

def printt(label):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(label)))

    return None

# 0. user defined variables
clean_fastq_dir='/home/adrian/projects/vigur/data/transcriptomics/tempo/'
threads = 20
results_dir='/home/adrian/projects/vigur/data/transcriptomics/tempo_kallisto/'
transcriptome_index='/home/adrian/software/kallisto/homo_sapiens/transcriptome.idx'
gtf_file = '/home/adrian/software/kallisto/homo_sapiens/Homo_sapiens.GRCh38.96.gtf'

# 1. recover labels
printt('recover labels...')

all_files = os.listdir(clean_fastq_dir)
labels = sorted(list(set(([element.split('1.R')[0] for element in all_files]))))
working_labels = labels[:3]

# 2. call kallisto quant
kallisto_caller(working_labels)
