import os
import sys
import argparse 
import numpy as np
from os import listdir, mkdir
from subprocess import call
from src.get_cluster import miltiprocess_analyze

def main():

    parser = argparse.ArgumentParser(description=' P I K E : metagenomic tool for noisy reads analyzing')
    group1 = parser.add_argument_group("main arguments")
    group2 = parser.add_argument_group("advance arguments")
    group3 = parser.add_argument_group("read quality arguments")
    group4 = parser.add_argument_group("clustering arguments")
    group5 = parser.add_argument_group("consensus building arguments")

    #______________________Main options_______________________________
    group1.add_argument('-fastq',
                        type=str,
                        help='Directory with merged read files',
                        required=True,
                        default=None)      
    group1.add_argument('-output',
                        type=str,
                        help='Directory with output files',
                        default='./Pike_results')
    #__________________________________________________________________
    #______________________Advance options_____________________________     
    group2.add_argument('-mode',
                        type=str,
                        help='Choice of analysis mode. single - each sample is analyzed separately, pool - joint analysis of all samples (single by default)',
                        required=True,
                        default='single')    
    group2.add_argument('-usereads', 
                        type=int,
                        help='Maximum number of reads per sample for analysis (Everything is used by default)',
                        default=np.inf)
    group2.add_argument('-threads', 
                        type=int,
                        help='Number of threads',
                        default=1)
    #__________________________________________________________________
    #______________________Read quality options________________________
    group3.add_argument('--trim_primer',
                        action='store_true',
                        help='Trim primers (If this mod is selected, then it is necessary to apply primer sequences)',
                        default=False)
    group3.add_argument('-primerF',
                        type=str,
                        help='Forward primer sequence',
                        default=None)   
    group3.add_argument('-primerR',
                        type=str,
                        help='Reverse primer sequence',
                        default=None)   
    group3.add_argument('-read_q_score',
                        type=int,
                        help='Median value of read quality',
                        default=20)
    group3.add_argument('-minlen',
                        type=str,
                        help='Minimum expected amplicon length (350 by default)',
                        default=350)
    group3.add_argument('-maxlen',
                        type=int,
                        help='Maximum expected amplicon length (600 by default)',
                        default=600)
    #_________________________________________________________________
    #______________________Clustering parametrs_______________________
    group4.add_argument('-umap_neighbours',
                        type=int,
                        help='Number of nearest neighbors for UMAP (30 by default)',
                        default=30)
    group4.add_argument('-cluster_size',
                        type=int,
                        help='Number of dots for HDBSCAN (30 by default)',
                        default=30)
    group4.add_argument('-k',
                        type=int,
                        help='k-mer size number (6 by default)',
                        default=6)
    group4.add_argument('--visualize', 
                        help='Create clustering figure (false by default)',
                        default=False,
                        action='store_true')
    #_________________________________________________________________
    #______________________Consensus parametrs________________________
    group5.add_argument('-consensus_seq_lim',
                        type=int,
                        help='Minimum number of reads supporting consensus sequence',
                        default=20)
    group5.add_argument('-letter_Q_lim',
                        type=int,
                        help='The minimum quality of the consensus letter. If the quality is less, then N is set (15 by default)',
                        default=15)
    #_________________________________________________________________

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    #Innitializing arguments 
    #____Main options__________________________
    path_to_fastq = args.fastq
    output = args.output
    #__________________________________________
    #____Advance options_______________________
    mode = args.mode
    usereads = args.usereads #Not for pool
    threads = args.threads
    #__________________________________________
    #____Read quality options__________________
    trim_primer = args.trim_primer
    primerF = args.primerF
    primerR = args.primerR
    minlen = args.minlen
    maxlen = args.maxlen
    read_q_score = args.read_q_score
    #__________________________________________
    #____Clustering parametrs__________________
    umap_neighbours = args.umap_neighbours
    hdbscan_neighbours = args.cluster_size
    k = args.k
    visualize = args.visualize
    #____Consensus parametrs__________________
    #__________________________________________
    consensus_seq_lim = args.consensus_seq_lim
    letter_Q_lim = args.letter_Q_lim
    #__________________________________________    


    try:        
        os.mkdir(output)
    except FileExistsError:

        print('The output directory already exists!')

    #output dir bases
    os.mkdir(f'{output}/read_preprocessing/')
    os.mkdir(f'{output}/work_dir/')
    os.mkdir(f'{output}/results/')
    #read processing dir dir bases
    if trim_primer == True:
        if primerF is None or primerR is None:
            print('For trimming, it is necessary to give primer sequences (primerF and primerR)!')
            sys.exit()

        os.mkdir(f'{output}/read_preprocessing/cutadapt_round1')
        os.mkdir(f'{output}/read_preprocessing/cutadapt_round2')

    os.mkdir(f'{output}/read_preprocessing/filtered_reads')    

    print('                      R U N   A N A L Y S I S                      ')
    print('===================================================================')

    miltiprocess_analyze(path_to_fastq, 
                         output,
                         mode,
                         threads,
                         usereads,
                         trim_primer,
                         primerF, 
                         primerR,
                         minlen, 
                         maxlen,
                         read_q_score,
                         umap_neighbours,
                         hdbscan_neighbours,
                         k,
                         visualize,
                         consensus_seq_lim,
                         letter_Q_lim)
    
    print('Your results are ready!')
    print('Thank you for using PIKE.')
    print('If you find a bug, please report it in a Git issues or contact by email.')

if __name__ == "__main__":
    main()