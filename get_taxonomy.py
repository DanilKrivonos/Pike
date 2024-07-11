import argparse 
import sys
from src.call_taxonomy import call_taxonomy
from pandas import read_csv

def main():

    parser = argparse.ArgumentParser(description=' P I K E taxonomy caller')

    #______________________Main options_______________________________
    parser.add_argument('-otutab',
                        type=str,
                        help='Path to otu table',
                        required=True,
                        default=None)      
    parser.add_argument('-dbpath',
                        type=str,
                        help='Path to taxonomy DB (SeqID Kingdom;Phylum;Class;Order;Family;Genus;Species)')
    
    parser.add_argument('-identity',
                        type=float,
                        help='Identity treshold (default : 95)',
                        default=95)
    parser.add_argument('-cov_lim',
                        type=float,
                        help='Caverage treshold (default : 60)',
                        default=60)
    parser.add_argument('-eval_lim',
                        type=float,
                        help='E-value treshold (default : 1e-05)',
                        default=60)
    parser.add_argument('-threads',
                        type=int,
                        help='Number of threads (default : 4)',
                        default=4)  
    parser.add_argument('-output',
                        type=str,
                        help='Directory with output files',
                        default='./Pike_results_tax')


    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    #Innitializing arguments 
    #__________________________Options__________________________
    otutab = args.otutab
    dbpath = args.dbpath
    identity = args.identity
    cov_lim = args.cov_lim
    eval_lim = args.eval_lim
    output = args.output
    threads = args.threads
    mearged_pike_out = read_csv(otutab, sep='\t', index_col=0)

    tax_table = call_taxonomy(output, 
                             dbpath,
                             mearged_pike_out.T, 
                             threads,
                             identity_filter=identity, 
                             cov_lim=cov_lim, 
                             evalue_filter=eval_lim)
    
    tax_table.to_csv(f'{output}/tax_table.tsv', sep='\t')

if __name__ == "__main__":
    main()