import numpy as np
import pandas as pd
from os import listdir
from pandas import read_csv, DataFrame

def merge_output(output, mode):

    merged_otu_table = []

    for samp in listdir(f'{output}/results/'):
        
        opn_res = read_csv(f'{output}/results/{samp}/results.tsv', sep='\t', index_col=0)
        count = 0
        merged_otu_table.append(DataFrame(data=opn_res.values, index=opn_res.index, columns=[samp]))
        
        if mode != 'pool':    
            with open(f'{output}/work_dir/{samp}/OTU.fasta', 'w') as opn_fasta:

                for line in opn_res.index:
                    
                    opn_fasta.write(f'>{samp}_{count}_{opn_res["Count"][line]}\n{line}\n')
                    count += 1
        else:
            with open(f'{output}/work_dir/{samp}_OTU.fasta', 'w') as opn_fasta:

                for line in opn_res.index:
                    
                    opn_fasta.write(f'>{samp}_{count}_{opn_res["Count"][line]}\n{line}\n')
                    count += 1

    merged_otu_table = pd.concat(merged_otu_table, axis=1).fillna(0)
    merged_otu_table = merged_otu_table[np.sort(merged_otu_table.columns)].T
    merged_otu_table.to_csv(f'{output}/merged_otu_table.tsv', sep='\t')
