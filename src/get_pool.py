import os
from os import listdir
import numpy as np
from umap import UMAP
from hdbscan import HDBSCAN
from subprocess import call
from sklearn.decomposition import PCA
from pandas import DataFrame, read_csv
from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from src.get_features import collect_features
from src.get_visualization import get_visualisation
from src.get_filtered_clusters import filter_cluster
from src.get_consensus import get_msa_info, get_consensus, medaka_run
from src.get_good_reads import run_trimming, run_filtering

def create_pool_dirs(output):

    #workdir dir bases
    os.mkdir(f'{output}/work_dir/features')
    os.mkdir(f'{output}/work_dir/proto_consensus')
    os.mkdir(f'{output}/work_dir/clusters_data_fastq')
    os.mkdir(f'{output}/work_dir/clusters_data_fasta')
    os.mkdir(f'{output}/work_dir/msa')
    os.mkdir(f'{output}/work_dir/medaka')

def collect_pool_features(interval,
                          path_to_fastq,
                          output,
                          usereads,
                          k,
                          trim_primer, 
                          primerF, 
                          primerR,
                          minlen, 
                          maxlen,
                          read_q_score):
                    
    for barcode in interval:
        
        barcode = barcode.replace('.fastq', '')
            
        print(' R E A D S   P R E P R O C E S S I N G ')
        print('=======================================')

        if trim_primer == True:
    
            path_to_fastq = run_trimming(path_to_fastq, 
                                         barcode, 
                                         output, 
                                         primerF, 
                                         primerR)

        path_to_fastq = run_filtering(path_to_fastq, 
                                      barcode, 
                                      output, 
                                      minlen, 
                                      maxlen, 
                                      read_q_score)

        print(' F E A T U R E S   C O L L E C T I N G ')
        print('=======================================')

        K_MERS_FREQ, GC_CONTENT, READ, READ_ID, LENS, QUALITY, BARCODE_ID = collect_features(path_to_fastq, 
                                                                                             barcode, 
                                                                                             usereads, 
                                                                                             k,
                                                                                             read_q_score) #Features collection
        save_tsv = DataFrame({'K_MERS_FREQ' : K_MERS_FREQ,
                              'GC_CONTENT' : GC_CONTENT, 
                              'READ' : READ,
                              'READ_ID' : READ_ID,
                              'LENS' : LENS, 
                              'QUALITY' : QUALITY, 
                              'BARCODE_ID' : BARCODE_ID
                              })

        save_tsv.to_csv(f'{output}/work_dir/features/{barcode}.tsv', sep='\t')
        
def run_pool(output,
             umap_neighbours,
             hdbscan_neighbours,
             visualize,
             consensus_seq_lim,
             letter_Q_lim):
    
    K_MERS_FREQ = [] 
    GC_CONTENT = []
    READ = []
    READ_ID =[] 
    LENS = []
    QUALITY = []
    BARCODE_ID = []
    #Barcodes features mearging
    for barcode_tsv in listdir(f'{output}/work_dir/features/'):

        opn_barcode_features = read_csv(f'{output}/work_dir/features/{barcode_tsv}', sep='\t', index_col=0)
        
        K_MERS_FREQ.extend([eval(read) for read in opn_barcode_features['K_MERS_FREQ'].values])
        GC_CONTENT.extend(list(opn_barcode_features['GC_CONTENT'].values))
        READ.extend([SeqRecord(read) for read in opn_barcode_features['READ'].values])
        READ_ID.extend(list(opn_barcode_features['READ_ID'].values))
        LENS.extend(list(opn_barcode_features['LENS'].values))
        QUALITY.extend(list(opn_barcode_features['QUALITY'].values))
        BARCODE_ID.extend(list(opn_barcode_features['BARCODE_ID'].values))

    #___Staged decomposition_________________________________________________________________________________________________________
        
    print('    C L U S T E R I N G   S T A G E    ')
    print('=======================================')

    pca_model = PCA(n_components=30,
                    random_state=0)
    pca_data = pca_model.fit_transform(K_MERS_FREQ)
    umap_model = UMAP(n_components=2,
                        n_neighbors=umap_neighbours,
                        min_dist=0.01,
                        metric='braycurtis',
                        random_state=0)
    umap_dat = umap_model.fit_transform(pca_data)
    #_______________________________________________________________________________________________________________________________
    
    RESULT_DICT = {'Read ID' : READ_ID,
                   'BARCODE' : BARCODE_ID, 
                   '1 UMAP COMPONENT' : umap_dat[:, 0], 
                   '2 UMAP COMPONENT' : umap_dat[:, 1], 
                   'Length' : LENS,
                   'fastq' : READ,
                   'GC content' : GC_CONTENT,
                   'QUALITY' : QUALITY,
                   'K-mers signature' : K_MERS_FREQ}
    RESULT_DF = DataFrame(RESULT_DICT) #All information mearged in df

    #___Cluster identification________________________________________________________________________________________________________
    hdbscan = HDBSCAN(min_cluster_size=hdbscan_neighbours,
                        cluster_selection_epsilon=0.01, 
                        gen_min_span_tree=True,
                        metric='braycurtis')
    clusters = hdbscan.fit_predict(umap_dat)
    RESULT_DF['Clusters'] = clusters
    #_______________________________________________________________________________________________________________________________

    RESULT_DF.to_csv(f'{output}/work_dir/metadata_df.tsv', sep='\t')
    clusters = RESULT_DF['Clusters'].unique()
    clusters = np.sort(clusters[clusters != -1])
    filtered_add = {}
    
    #___Cluster processing_____________________________________________________________________________________________________________
    for clt in clusters:
        if clt == -1:
            continue
        clt_dat = RESULT_DF[RESULT_DF['Clusters'] == clt]
        filtered_data = filter_cluster(clt_dat) #filtering clustaer data
        
        if len(filtered_data) < consensus_seq_lim: # can be 0 length for artefact clusters
            continue

        filtered_add[clt] = filtered_data

        #_______Recording cluster data_________________________________________________________________________________________________
                
        clt_fastq = open(f'{output}/work_dir/clusters_data_fastq/{clt}.fastq', 'w')
        clt_fasta = open(f'{output}/work_dir/clusters_data_fasta/{clt}.fasta', 'w')
        quality_dict = {}
        
        for idx in filtered_data.index:

            read_record = filtered_data.fastq[idx]
            clt_fastq.write(read_record.format('fastq')) #fastq recording
            quality_dict[read_record.id] = [i for i in read_record.format('fastq').split('\n')[3]] #letter quality adding
            clt_fasta.write(f'>{read_record.id}\n{read_record.seq}\n') #fasta recording
        
        clt_fastq.close()
        clt_fasta.close()
        #_______________________________________________________________________________________________________________________________
        
        #_______Cluster Consensus Building______________________________________________________________________________________________
                
        call(f'mafft --quiet {output}/work_dir/clusters_data_fasta/{clt}.fasta > {output}/work_dir/msa/{clt}.msa', shell=True)
        consensus_possitions, possitional_Q_distributions = get_msa_info(output, 
                                                                         barcode, 
                                                                         clt, 
                                                                         quality_dict) #MSA info collection
        protoconsensus = get_consensus(consensus_possitions, 
                                        possitional_Q_distributions, 
                                        letter_Q_lim) #Building a protoconsensus

        with open(f'{output}/work_dir/proto_consensus/{clt}.fasta', 'w') as protocons:

            protocons.write(f'>clt_{clt}_{len(clt_dat)}\n{protoconsensus}\n') #Recording a protoconsensus    
        #_______________________________________________________________________________________________________________________________

        #_______Building visualization of clustering___________________________________________________________________________________
        
        if visualize == True:
            
            get_visualisation(output, barcode, filtered_add)
    #_______________________________________________________________________________________________________________________________
    
    print('    P O L I S H   C O N S E N S U S    ')
    print('=======================================')
    
    medaka_run(output, '/')
    #demultiplex by barcodes
    demultiplex_results(output)

def demultiplex_results(output):
        
    print('  O U T P U T   P R E P A R A T I O N  ')
    print('=======================================')

    otu_decoder = {}

    for clt in listdir(f'{output}/work_dir/medaka/'):

        opn_consensus = parse(f'{output}/work_dir/medaka/{clt}/consensus.fasta', 'fasta')

        for line in opn_consensus:
            
            otu_decoder[clt] = line.seq
    
    metadata_df = read_csv(f'{output}/work_dir/metadata_df.tsv', sep='\t')

    for barcode in metadata_df['BARCODE'].unique():

        os.mkdir(f'{output}/results/{barcode}')

        barcode_df = metadata_df[metadata_df['BARCODE'] == barcode]
        representation = barcode_df['Clusters'].value_counts('Length').to_dict()

        barcode_results = DataFrame({otu_decoder[clt] : representation[clt] for clt in representation.keys()})
        barcode_results.to_csv(f'{output}/results/{barcode}/results.tsv', sep='\t')

