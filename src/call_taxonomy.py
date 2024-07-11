import os
from Bio.SeqIO import parse
from subprocess import call
from pandas import read_csv, DataFrame


def create_fasta(output, mearged_pike_out, dbpath):
    # Create fasta file 
    consensus = {}
    cons_conter = 0
    
    with open(f'{output}/all_consensus.fasta', 'w') as opn_fasta:
        for cons in mearged_pike_out.index:
    
            opn_fasta.write(f'>{cons_conter}\n{cons}\n')
            consensus[cons_conter] = cons
            cons_conter += 1
    
    return consensus
    
def run_blast(base, path, threads):
    
    call(f'makeblastdb -in {base} -dbtype nucl', shell=True)
    call(f'blastn -num_threads {threads}  -outfmt "7 qseqid sseqid pident evalue qcovs bitscore" -query {path}/all_consensus.fasta  -db {base} -out {path}/blast_results.txt', shell=True)

def decode_tax(base) -> dict:
    
    # DB decoder 
    # Use db header format: Kingdom    Phylum    Class    Order    Family    Genus    Species
    
    base = parse(base, 'fasta')
    taxonomy_linage = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    tax_decoder = {}
    
    for line in base:
        
        tax_decoder[line.id] = {}
        linage = line.description.split(';')
        linage[0] = linage[0].split()[1]
    
        for i in range(len(taxonomy_linage)):
            try:
                
                tax_decoder[line.id][taxonomy_linage[i]] = linage[i]
            
            except:
                
                tax_decoder[line.id][taxonomy_linage[i]] = 'NA'
    
    return tax_decoder

def parse_blast(path, 
                base, 
                data_tax, 
                consensus, 
                identity_filter, 
                cov_lim, 
                evalue_filter):
    
    # parser of blast table
    
    blast_header = ['qseqid',
                    'sseqid', 
                    'pident',
                    'evalue',
                    'qcovs', 
                    'bitscore']
    
    blasting_results = {}
    opn_blast = read_csv(f'{path}/blast_results.txt', sep='\t', comment='#', header=None, names=blast_header)
    
    for i in opn_blast['qseqid'].unique():
        
        blast_subset = opn_blast[opn_blast["qseqid"] == i]
        blast_subset = blast_subset[blast_subset['pident'] >= identity_filter]
        blast_subset = blast_subset[blast_subset['evalue'] <= evalue_filter]
        blast_subset = blast_subset[blast_subset['qcovs'] >= cov_lim]

        blast_subset = blast_subset.sort_values(by='evalue')
        blast_subset = blast_subset.sort_values(by='pident')[::-1]

        if len(blast_subset['sseqid'].values) == 0:
            continue
            
        subject = blast_subset['sseqid'].values[0]
        blasting_results[consensus[i]] = data_tax[subject]
        
    blasting_results_df = DataFrame(blasting_results).T
    
    return blasting_results_df


def call_taxonomy(output, 
                dbpath,
                mearged_pike_out, 
                threads,
                identity_filter=95, 
                cov_lim=60, 
                evalue_filter=1e-05):

    # Creating output directory
    try:        
        os.mkdir(output)
    except FileExistsError:
        print('The output directory already exists!')
        
    consensus = create_fasta(output, mearged_pike_out, dbpath)
    run_blast(dbpath, output, threads)
    data_tax = decode_tax(dbpath)
    blasting_results_df = parse_blast(output, 
                                    dbpath, 
                                    data_tax, 
                                    consensus, 
                                    identity_filter, 
                                    cov_lim, 
                                    evalue_filter)

    return blasting_results_df

