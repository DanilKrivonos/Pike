import os
from subprocess import call

def run_trimming(path_to_fastq, 
                 sample, 
                 output, 
                 primerF, 
                 primerR):

    command1 = f'cutadapt -g {primerF} --rc --discard-untrimmed --quiet -e 30 -o {output}/read_processing/cutadapt_round1/{sample}.fastq {path_to_fastq}/{sample}.fastq'
    command2 = f'cutadapt -g {primerR} --rc --discard-untrimmed --quiet -e 30 -o {output}/read_processing/cutadapt_round2/{sample}.fastq {output}/read_processing/cutadapt_round1/{sample}.fastq'
    #Trimming
    call(command1, shell=True)
    call(command2, shell=True)
    goodreads = f'{output}/read_processing/cutadapt_round2/'

    return goodreads

def run_filtering(path_to_fastq, 
                  sample, 
                  output, 
                  minlen, 
                  maxlen, 
                  read_q_score):

    command = f'filtlong --min_length {minlen} --max_length {maxlen} --mean_q_weight {read_q_score} {path_to_fastq}/{sample}.fastq > {output}/read_processing/filtered_reads/{sample}.fastq'
    #Filtering
    call(command, shell=True)
    goodreads = f'{output}/read_processing/filtered_reads/'

    return goodreads