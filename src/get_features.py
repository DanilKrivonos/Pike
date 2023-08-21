import numpy as np
from itertools import product 
from Bio.SeqIO import parse

def compress_read(cutted_seq):
    
    new_seq = ''

    for idx in range(len(cutted_seq) - 1):
        if cutted_seq[idx] == cutted_seq[idx + 1]:
            continue

        new_seq += cutted_seq[idx]

    if cutted_seq[idx + 1] != new_seq[-1]:

        new_seq += cutted_seq[idx + 1]
        
    return new_seq 


def normalize(kmers):
    
    norm = sum(list(kmers.values()))
    
    for kmer in kmers.keys():
        
        kmers[kmer] = kmers[kmer]/ norm
    
    return kmers   

def get_kmers_signature(seq, k=6):
    
    kmers = {"".join(kmer) : 0 for kmer in list(product("AGTC", repeat=k))}
    step = 1
    start = 0
    end = k
    reverse_seq = compress_read(seq.reverse_complement())
    foraverd_seq = compress_read(seq)
    
    while end != len(foraverd_seq) - 1:
        
        subseq = seq[start: end]         
        kmers[subseq] += 1
        
        revers_subseq = reverse_seq[start: end]
        kmers[revers_subseq] += 1
        
        start, end = start + step, end + step

    kmers = normalize(kmers)
    
    return kmers

def collect_features(path_to_fastq, 
                     barcode, 
                     usereads, 
                     k,
                     median_Q_lim):

    open_fastq = parse(f'{path_to_fastq}/{barcode}.fastq' , 'fastq')

    K_MERS_FREQ = []
    GC_CONTENT = []
    READ_ID = []
    READ_Seq = []
    READ_Q = []
    LENS = []
    QUALITY = []
    BARCODE_ID = [] #For pool mode
    read_counter = 0
    Q_have = 0
    
    #___Features collection______________________________________________________________
    for seq in open_fastq:
        read_counter += 1
        if Q_have == usereads:
            break
        quallist = np.array(list(map(ord,list(seq.format('fastq').split('\n')[-2]))))-33
        quallist = quallist[quallist < 90]
        median_Q_score = np.median(quallist)
        
        if median_Q_score < median_Q_lim:
            continue
            
        GC_count = str(seq.seq).count('G') + str(seq.seq).count('C')
        K_MERS_FREQ.append(list(get_kmers_signature(seq.seq, k=k).values()))
        GC_CONTENT.append(GC_count/ len(seq.seq))
        READ_ID.append(seq.id)
        READ_Seq.append(seq.seq)
        READ_Q.append(seq.format('fastq').split('\n')[3])
        LENS.append(len(seq.seq))
        QUALITY.append(median_Q_score)
        BARCODE_ID.append(barcode)
        Q_have += 1
    #_____________________________________________________________________________________
    print(f'{np.round((1-Q_have/read_counter)*100)}% or reads was dropped') #Добавить log
    print(f'Will be used {len(K_MERS_FREQ)}')
    
    return K_MERS_FREQ, GC_CONTENT, READ_ID, READ_Seq, READ_Q, LENS, QUALITY, BARCODE_ID

        