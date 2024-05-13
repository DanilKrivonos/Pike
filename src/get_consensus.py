import os
import numpy as np
from os import listdir
from pandas import DataFrame
from Bio.SeqIO import parse
from subprocess import call

def get_msa_info(output, 
                 barcode, 
                 clt,
                 quality_dict):
    
    opn_msa = parse(f'{output}/work_dir/{barcode}/msa/{clt}.msa', 'fasta')
    consensus_possitions = {}
    possitional_Q_distributions = {}

    for line in opn_msa:
        
        shift = 0
        
        for possition in range(len(line.seq)):
            
            if possition not in consensus_possitions:
    
                consensus_possitions[possition] = []
                possitional_Q_distributions[possition] = {}
            consensus_possitions[possition].append(line.seq[possition])
            
            if line.seq[possition] == '-':
                shift += 1
                continue
            if line.seq[possition] not in possitional_Q_distributions[possition]:
                
                possitional_Q_distributions[possition][line.seq[possition]] = []
                
            possitional_Q_distributions[possition][line.seq[possition]].append(ord(quality_dict[line.id][possition - shift])-33)
    
    return consensus_possitions, possitional_Q_distributions

def most_common(lst):

    return max(set(lst), key=lst.count)

def get_consensus(consensus_possitions, 
                  possitional_Q_distributions, 
                  letter_Q_lim=15):
    
    protoconsensus = ""
    
    for i in consensus_possitions.keys():

        letter_possition = most_common(consensus_possitions[i])
        
        if letter_possition != '-':
        
            max_qaul = 0
            letter_possition = 'N'
            
            for alt_pos in possitional_Q_distributions[i].keys():
                
                if len(possitional_Q_distributions[i][alt_pos]) <= 3:
                    continue
                alt_qual = np.mean(possitional_Q_distributions[i][alt_pos])#/np.std(qal_dist[i][alt_pos])
            
                if alt_qual > max_qaul:
                    
                    letter_possition = alt_pos
                    max_qaul = alt_qual
                    
            if max_qaul < letter_Q_lim:
    
                letter_possition = 'N'
                
        protoconsensus += letter_possition
     
    protoconsensus = protoconsensus.replace('-', '').upper()
    
    return protoconsensus


def medaka_run(output, barcode, threads=4):
    
    for cluster in listdir(f'{output}/work_dir/{barcode}/clusters_data_fastq/'):
        
        cluster = cluster.split('.')[0]
        reads = f'{output}/work_dir/{barcode}/clusters_data_fastq/{cluster}.fastq'
        protoconsensus = f'{output}/work_dir/{barcode}/proto_consensus/{cluster}.fasta'
        out_dir = f'{output}/work_dir/{barcode}/medaka/{cluster}'
        command = f'medaka_consensus -t {threads} -i {reads} -d {protoconsensus} -o {out_dir}'
        call(command, shell=True)

def trim_N(seq, limit=0.1):
    
    cut_start = 0
    cut_end = len(seq)
    
    for idx in range(len(seq)):
        if seq[idx] == 'N':
            if idx <= len(seq) * limit:
            
                cut_start = idx
    
        if seq[idx] == 'N':
            if idx >= len(seq)*(1-limit):
                
                cut_end = idx
                break
    
    trimmed_seq = seq[cut_start + 1: cut_end]
    
    return trimmed_seq

def prepare_output(output, barcode):
    
    otu_dict = {}
    
    for clt in listdir(f'{output}/work_dir/{barcode}/medaka'):
        
        opn_consensus = parse(f'{output}/work_dir/{barcode}/medaka/{clt}/consensus.fasta', 'fasta')
        
        for line in opn_consensus:
            
            otu = str(line.seq)
            otu = trim_N(otu)
            otu_count = int(line.id.split('_')[-1])
            
            if otu not in otu_dict.keys():
            
                otu_dict[otu] = {'Count' : 0}
                
            otu_dict[otu]['Count'] += otu_count
    
    otu_df = DataFrame(otu_dict).T
    otu_df.to_csv(f'{output}/results/{barcode}/results.tsv', sep='\t')