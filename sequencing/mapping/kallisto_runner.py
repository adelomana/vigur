###
### usage: time python kallisto_runner.py &> messages.txt
###

import sys, datetime, os

def kallisto_caller(label):

    printt('about to quantify {}'.format(label))

    sample_output_dir = results_dir + label
    executable = 'time kallisto quant'
    options = ' -i {} -o {} --bias -t {} -b {} {} '.format(transcriptome_index, sample_output_dir, threads, boots, strand_flag)
    fastq_files = '{} {}'.format(clean_fastq_dir + label + '/' + label + '_R1_clean.fastq.gz', clean_fastq_dir + label + '/' + label + '_R2_clean.fastq.gz')
    command = executable + options + fastq_files

    print('')
    print(command)
    os.system(command)
    print('')

    return None

def printt(label):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(label)))

    return None

###
### 0. user defined variables
###
clean_fastq_dir = '/home/adrian/scratch/decode2/clean_fastq/'
boots = 100
threads = 20
results_dir = '/home/adrian/projects/vigur/results/sequencing/kallisto/kallisto.{}/'.format(boots)
transcriptome_index = '/home/adrian/software/kallisto/ensembl_v96/transcriptome.idx'

strand_flag = '--rf-stranded'   # processed 115,770,977 reads, 51,027,094 reads pseudoaligned
strand_flag = '--fr-stranded'   # processed 115,770,977 reads, 50,260,575 reads pseudoaligned
strand_flag = ''                # processed 115,770,977 reads, 101,434,048 reads pseudoaligned

###
### 1. recover labels
###
printt('recover labels...')

labels = os.listdir(clean_fastq_dir)
labels.sort()

###
### 2. call kallisto quant
###
if os.path.exists(results_dir) == False:
    os.mkdir(results_dir)

for label in labels:
    kallisto_caller(label)
