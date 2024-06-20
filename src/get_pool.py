import os
import pickle
from os import listdir
import numpy as np
from umap import UMAP
from umap.parametric_umap import ParametricUMAP
from hdbscan import HDBSCAN
from subprocess import call
from skbio.stats.composition import clr
from sklearn.decomposition import PCA
from pandas import DataFrame, read_csv
from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from src.get_features import collect_features
from src.get_filtered_clusters import filter_cluster
from src.get_consensus import get_msa_info, get_consensus, medaka_run
from src.get_good_reads import run_trimming, run_filtering
from scipy.sparse import vstack, csr_matrix

def create_pool_dirs(output):

    #workdir dir bases
    os.mkdir(f'{output}/work_dir/features')
    os.mkdir(f'{output}/work_dir/features/tsv')
    os.mkdir(f'{output}/work_dir/features/sparse')
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
        
        barcode_path = path_to_fastq

        barcode = barcode.replace('.fastq', '')
            
        print(' R E A D S   P R E P R O C E S S I N G ')
        print('=======================================')

        if trim_primer == True:
    
            barcode_path = run_trimming(barcode_path, 
                                         barcode, 
                                         output, 
                                         primerF, 
                                         primerR)

        barcode_path = run_filtering(barcode_path, 
                                      barcode, 
                                      output, 
                                      minlen, 
                                      maxlen, 
                                      read_q_score)

        print(' F E A T U R E S   C O L L E C T I N G ')
        print('=======================================')

        K_MERS_FREQ, GC_CONTENT, READ_ID, READ_Seq, READ_Q, LENS, QUALITY, BARCODE_ID = collect_features(barcode_path, 
                                                                                                         barcode, 
                                                                                                         usereads, 
                                                                                                         k,
                                                                                                         read_q_score) #Features collection
        save_tsv = DataFrame({'Read ID' : READ_ID,
                              'BARCODE' : BARCODE_ID, 
                              'Length' : LENS,
                              'READ_Q' : READ_Q,
                              'READ_Seq' : READ_Seq, 
                              'GC content' : GC_CONTENT,
                              'QUALITY' : QUALITY})#,
                              #'K-mers signature' : list(K_MERS_FREQ)})
        
        with open(f'{output}/work_dir/features/sparse/K_MERS_FREQ_{barcode}.obj', 'wb') as fp:

            pickle.dump(K_MERS_FREQ, fp)

        save_tsv.to_csv(f'{output}/work_dir/features/tsv/{barcode}.tsv', sep='\t')

def demultiplex_results(output):
        
    print('  O U T P U T   P R E P A R A T I O N  ')
    print('=======================================')

    otu_decoder = {}
    good_clt_list = []

    for clt in listdir(f'{output}/work_dir/medaka/'):

        good_clt_list.append(int(clt))
        opn_consensus = parse(f'{output}/work_dir/medaka/{clt}/consensus.fasta', 'fasta')

        for line in opn_consensus:
            
            otu_decoder[int(clt)] = line.seq
    
    metadata_df = read_csv(f'{output}/work_dir/metadata_df.tsv', sep='\t', index_col=0)
    metadata_df = metadata_df.loc[metadata_df['Clusters'].isin(good_clt_list)]

    for barcode in metadata_df['BARCODE'].unique():

        os.mkdir(f'{output}/results/{barcode}')

        barcode_df = metadata_df[metadata_df['BARCODE'] == barcode]
        representation = barcode_df.value_counts('Clusters').to_dict()
        barcode_results_dict = {otu_decoder[clt] : {'Count' : representation[clt]} for clt in representation.keys()}
        barcode_results = DataFrame(barcode_results_dict).T
        barcode_results.to_csv(f'{output}/results/{barcode}/results.tsv', sep='\t')


def run_pool(output,
             threads,
             umap_neighbours,
             hdbscan_neighbours,
             consensus_seq_lim,
             letter_Q_lim):
    
    K_MERS_FREQ = [] 
    GC_CONTENT = []
    READ_Q = []
    READ_Seq = []
    READ_ID =[] 
    LENS = []
    QUALITY = []
    BARCODE_ID = []
    #Barcodes features mearging
    for barcode_tsv in listdir(f'{output}/work_dir/features/tsv'):

        opn_barcode_features = read_csv(f'{output}/work_dir/features/tsv/{barcode_tsv}', sep='\t', index_col=0)
        barcode = barcode_tsv.split('.')[0]

        with open(f'{output}/work_dir/features/sparse/K_MERS_FREQ_{barcode}.obj', 'rb') as fp:
            
            K_MERS_FREQ_barcode = pickle.load(fp)
        

        K_MERS_FREQ.append(K_MERS_FREQ_barcode)
        GC_CONTENT.extend(list(opn_barcode_features['GC content'].values))
        READ_Q.extend(list(opn_barcode_features['READ_Q'].values))
        READ_Seq.extend(list(opn_barcode_features['READ_Seq'].values))
        READ_ID.extend(list(opn_barcode_features['Read ID'].values))
        LENS.extend(list(opn_barcode_features['Length'].values))
        QUALITY.extend(list(opn_barcode_features['QUALITY'].values))
        BARCODE_ID.extend(list(opn_barcode_features['BARCODE'].values))

    K_MERS_FREQ = vstack(K_MERS_FREQ)
    print(K_MERS_FREQ.A.shape)
    #___Staged decomposition_________________________________________________________________________________________________________
        
    print('    C L U S T E R I N G   S T A G E    ')
    print('=======================================')
    #K_MERS_FREQ = np.array(K_MERS_FREQ)
    #idx = np.argwhere(np.all(K_MERS_FREQ[..., :] == 0, axis=0))
    #K_MERS_FREQ = np.delete(K_MERS_FREQ, idx, axis=1)
    #K_MERS_FREQ += 1
    #normalize = lambda x: x / np.sum(x)
    #K_MERS_FREQ = np.array(list(map(normalize, K_MERS_FREQ)))
    clr_data = csr_matrix(clr(K_MERS_FREQ.A))#clr_data = clr(K_MERS_FREQ)
    
    pca_model = PCA(n_components=15,
                    random_state=0, svd_solver='arpack')
    pca_data = pca_model.fit_transform(clr_data)
    # umap_model = UMAP(n_components=2,
    #                  n_neighbors=umap_neighbours,
    #                  min_dist=0.01,
    #                  metric='euclidean',
    #                  random_state=0)
    umap_model = ParametricUMAP(n_components=2)
    umap_dat = umap_model.fit_transform(pca_data)
    #_______________________________________________________________________________________________________________________________
    
    RESULT_DICT = {'Read ID' : READ_ID,
                    'BARCODE' : BARCODE_ID, 
                    '1 UMAP COMPONENT' : umap_dat[:, 0], 
                    '2 UMAP COMPONENT' : umap_dat[:, 1], 
                    'Length' : LENS,
                    'READ_Q' : READ_Q,
                    'READ_Seq' : READ_Seq, 
                    'GC content' : GC_CONTENT,
                    'QUALITY' : QUALITY}#,
                    #'K-mers signature' : list(K_MERS_FREQ)}
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

            read_record = f'@{filtered_data["Read ID"][idx]}\n{filtered_data["READ_Seq"][idx]}\n+\n{filtered_data["READ_Q"][idx]}\n'
            fasta_record = f'>{filtered_data["Read ID"][idx]}\n{filtered_data["READ_Seq"][idx]}\n'
            clt_fastq.write(read_record) #fastq recording
            quality_dict[filtered_data["Read ID"][idx]] = [i for i in RESULT_DF['READ_Q'][idx]] #letter quality adding
            clt_fasta.write(fasta_record) #fasta recording
        
        clt_fastq.close()
        clt_fasta.close()
        #_______________________________________________________________________________________________________________________________
        
        #_______Cluster Consensus Building______________________________________________________________________________________________
                
        call(f'mafft --thread {threads} --quiet {output}/work_dir/clusters_data_fasta/{clt}.fasta > {output}/work_dir/msa/{clt}.msa', shell=True)
        consensus_possitions, possitional_Q_distributions = get_msa_info(output, 
                                                                         '/', 
                                                                         clt, 
                                                                         quality_dict) #MSA info collection
        protoconsensus = get_consensus(consensus_possitions, 
                                        possitional_Q_distributions, 
                                        letter_Q_lim) #Building a protoconsensus

        with open(f'{output}/work_dir/proto_consensus/{clt}.fasta', 'w') as protocons:

            protocons.write(f'>clt_{clt}_{len(clt_dat)}\n{protoconsensus}\n') #Recording a protoconsensus    
        #_______________________________________________________________________________________________________________________________

    #_______________________________________________________________________________________________________________________________
    
    print('    P O L I S H   C O N S E N S U S    ')
    print('=======================================')
    
    medaka_run(output, '/', threads)
    #demultiplex by barcodes
    demultiplex_results(output)