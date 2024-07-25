import numpy as np
from itertools import product 
from Bio.SeqIO import parse
from scipy.sparse import csr_array, vstack


def compress_read(cutted_seq):
    
    new_seq = ''

    for idx in range(len(cutted_seq) - 1):
        if cutted_seq[idx] == cutted_seq[idx + 1]:
            continue

        new_seq += cutted_seq[idx]

    if cutted_seq[idx + 1] != new_seq[-1]:

        new_seq += cutted_seq[idx + 1]
        
    return new_seq 

def get_kmers_signature(seq, k=6):
    
    kmers = {"".join(kmer) : 0 for kmer in list(product("AGTCN", repeat=k))}
 #   step = 1
  #  start = 0
  #  end = k
    reverse_seq = compress_read(seq.reverse_complement())
    forward_seq = compress_read(seq)
    
    #limit = len(forward_seq) - 1

    for start in  range(len(forward_seq) - k):

        end = start + k
        subseq = forward_seq[start: end]
        kmers[subseq] += 1
    #    revers_subseq = reverse_seq[start: end]
    #    kmers[revers_subseq] += 1
    #while end != limit:
        
     #   subseq = forward_seq[start: end]         
     #   kmers[subseq] += 1
        
     #   revers_subseq = reverse_seq[start: end]
     #   kmers[revers_subseq] += 1
        
      #  start, end = start + step, end + step
    
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
        
    #    kmer_res = np.array(list(get_kmers_signature(seq.seq, k=k).values()))
      #  kmer_res += 1
      #  kmer_res = kmer_res / np.sum(kmer_res)

#        K_MERS_FREQ.append(csr_array(kmer_res))
        kmer_res = list(get_kmers_signature(seq.seq, k=k).values())

        K_MERS_FREQ.append(kmer_res)
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
 #   K_MERS_FREQ = vstack(K_MERS_FREQ)

    return K_MERS_FREQ, GC_CONTENT, READ_ID, READ_Seq, READ_Q, LENS, QUALITY, BARCODE_ID

        