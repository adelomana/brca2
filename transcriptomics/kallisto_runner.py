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

boots = 100
threads = 20

clean_fastq_dir = '/home/adrian/projects/brca2/results/clean_fastq/'
results_dir = '/home/adrian/projects/brca2/results/kallisto/kallisto.{}/'.format(boots)
transcriptome_index = '/home/adrian/software/kallisto/ensembl_v96/transcriptome.idx'

#strand_flag = '--rf-stranded'
#strand_flag = '--fr-stranded'
strand_flag = ''

###
### 1. recover labels
###
printt('recover labels...')

labels = os.listdir(clean_fastq_dir)
labels.sort()
print(labels)

###
### 2. call kallisto quant
###
if os.path.exists(results_dir) == False:
    os.mkdir(results_dir)

for label in labels:
    kallisto_caller(label)
    #sys.exit()
