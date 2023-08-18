import os
from subprocess import call

def run_trimming(path_to_fastq, 
                 sample, 
                 output, 
                 primerF, 
                 primerR):
    
    #Trimming round 1
    command_input1 = f'{path_to_fastq}/{sample}.fastq'
    command_output1 = f'{output}/read_preprocessing/cutadapt_round1/{sample}.fastq'    
    command1 = f'cutadapt -g {primerF} --rc --discard-untrimmed --quiet -e 30 -o {command_output1} {command_input1}'
    call(command1, shell=True)
    #Trimming round 2
    command_input2 = f'{output}/read_preprocessing/cutadapt_round1/{sample}.fastq'
    command_output2 = f'{output}/read_preprocessing/cutadapt_round2/{sample}.fastq'
    command2 = f'cutadapt -g {primerR} --rc --discard-untrimmed --quiet -e 30 -o {command_output2} {command_input2}'
    call(command2, shell=True)
    #new path to preprocessed fastq
    goodreads = f'{output}/read_preprocessing/cutadapt_round2/'

    return goodreads

def run_filtering(path_to_fastq, 
                  sample, 
                  output, 
                  minlen, 
                  maxlen, 
                  read_q_score):
    
    #Filtering
    command_input = f'{read_q_score} {path_to_fastq}/{sample}.fastq'
    command_output = f'{output}/read_preprocessing/filtered_reads/{sample}.fastq'
    command = f'filtlong --min_length {minlen} --max_length {maxlen} --mean_q_weight {command_input} > {command_output}'
    call(command, shell=True)
    #new path to preprocessed fastq
    goodreads = f'{output}/read_preprocessing/filtered_reads/'

    return goodreads