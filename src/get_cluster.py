import os
import random
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from os import listdir
from umap import UMAP
from hdbscan import HDBSCAN
from subprocess import call
from pandas import DataFrame
from multiprocessing import Process
from sklearn.decomposition import PCA
from src.get_pool import run_pool
from src.get_features import collect_features
from scipy.spatial.distance import pdist, squareform
from src.get_good_reads import run_trimming, run_filtering
from src.get_consensus import get_msa_info, get_consensus, medaka_run, prepare_output

def create_dirs(output, 
                barcode):

    #workdir dir bases
    os.mkdir(f'{output}/work_dir/{barcode}')
    os.mkdir(f'{output}/results/{barcode}')
    os.mkdir(f'{output}/work_dir/{barcode}/proto_consensus')
    os.mkdir(f'{output}/work_dir/{barcode}/clusters_data_fastq')
    os.mkdir(f'{output}/work_dir/{barcode}/clusters_data_fasta')
    os.mkdir(f'{output}/work_dir/{barcode}/msa')
    os.mkdir(f'{output}/work_dir/{barcode}/medaka')


def filter_cluster(clt_dat):
    
    #_____________FILTERING_BY_LENGTH_____________________________________
    median = np.median(clt_dat['Length'].values)
    st_div = np.std(clt_dat['Length'].values)
    left = median - 10
    right = median + 10
    filtered_data = clt_dat[clt_dat['Length'] > left]
    filtered_data = filtered_data[filtered_data['Length'] < right]
    #_____________________________________________________________________
    #_______________FILTERING_BY_GC_______________________________________
    median = np.median(clt_dat['GC content'].values)
    st_div = np.std(clt_dat['GC content'].values)
    left = median - 2*st_div
    right = median + 2*st_div
    filtered_data = filtered_data[filtered_data['GC content'] > left]
    filtered_data = filtered_data[filtered_data['GC content'] < right]
    #_____________________________________________________________________
    
    return  filtered_data

def get_color(tax_list):
    
    color = ''
    
    while color not in tax_list.values() and color == '':
        
        color = "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
    
    return color

def get_visualisation(output, 
                      barcode, 
                      filtered_add):

    fig, axs = plt.subplots(figsize=(10, 10))
    colordict = {}
    
    for clt in filtered_add.keys():
        
        mean_dists = np.median(squareform(pdist(filtered_add[clt][['1 UMAP COMPONENT', '2 UMAP COMPONENT']].values), 'braucyrtis'), axis=1)
        ref_idx = filtered_add[clt][mean_dists == np.min(mean_dists)].index[0]

        #____Looking for cluster centroids for figure______
        x =  filtered_add[clt].loc[ref_idx]['1 UMAP COMPONENT']
        y = filtered_add[clt].loc[ref_idx]['2 UMAP COMPONENT']
        colordict[clt] = get_color(colordict)
        
        plt.text(x, y, 
                 f'Cluster {clt}', 
                 fontsize=3, 
                 fontweight='book')
        sns.scatterplot(data=filtered_add[clt], 
                        x='1 UMAP COMPONENT', 
                        y='2 UMAP COMPONENT', 
                        s=1,
                        alpha=.3,
                        color=colordict[clt],
                        ax=axs)
        sns.scatterplot(data=filtered_add[clt], 
                        x='1 UMAP COMPONENT', 
                        y='2 UMAP COMPONENT', 
                        s=3,
                        alpha=.7,
                        color=colordict[clt],
                        ax=axs)
    #sns.despine(offset=10, trim=True)
    plt.xlabel('1 UMAP COMPONENT', fontweight='book')
    plt.ylabel('2 UMAP COMPONENT', fontweight='book')
    plt.grid()
    plt.legend([],[], frameon=False)
    plt.savefig(f'{output}/{barcode}/reults/{barcode}.pdf')
    plt.savefig(f'{output}/{barcode}/reults/{barcode}.pdf', dpi=800)

def run_barcodes(interval,
                 path_to_fastq,
                 output,
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
                 letter_Q_lim):
                    
    for barcode in interval:
        
        barcode = barcode.replace('.fastq', '')
        
        create_dirs(output, barcode) #creating and removing directoryes
        
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

        RESULT_DF.to_csv(f'{output}/work_dir/{barcode}/metadata_df.tsv', sep='\t')
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
                    
            clt_fastq = open(f'{output}/work_dir/{barcode}/clusters_data_fastq/{clt}.fastq', 'w')
            clt_fasta = open(f'{output}/work_dir/{barcode}/clusters_data_fasta/{clt}.fasta', 'w')
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
                    
            call(f'mafft --quiet {output}/work_dir/{barcode}/clusters_data_fasta/{clt}.fasta > {output}/work_dir/{barcode}/msa/{clt}.msa', shell=True)
            consensus_possitions, possitional_Q_distributions = get_msa_info(output, 
                                                                             barcode, 
                                                                             clt, 
                                                                             quality_dict) #MSA info collection
            protoconsensus = get_consensus(consensus_possitions, 
                                           possitional_Q_distributions, 
                                           letter_Q_lim) #Building a protoconsensus
    
            with open(f'{output}/work_dir/{barcode}/proto_consensus/{clt}.fasta', 'w') as protocons:
    
                protocons.write(f'>clt_{clt}_{len(clt_dat)}\n{protoconsensus}\n') #Recording a protoconsensus    
            #_______________________________________________________________________________________________________________________________

            #_______Building visualization of clustering___________________________________________________________________________________
            
            if visualize == True:
                
                get_visualisation(output, barcode, filtered_add)
        #_______________________________________________________________________________________________________________________________
        
        print('    P O L I S H   C O N S E N S U S    ')
        print('=======================================')
        
        medaka_run(output, barcode)

        print('  O U T P U T   P R E P A R A T I O N  ')
        print('=======================================')

        prepare_output(output, barcode)

def miltiprocess_analyze(path_to_fastq, 
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
                         letter_Q_lim):
        
    barcodes = listdir(path_to_fastq)
    intervals = np.linspace(0, len(barcodes)+1, threads+1)
    intervals = [barcodes[int(intervals[i]):int(intervals[i+1])] for i in range(len(intervals)-1)] # split to processes 
    procs = []

    if mode == 'single':
        for interval in intervals:
            
            proc = Process(target=run_barcodes, 
                        args=(interval,
                              path_to_fastq,
                              output,
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
                              letter_Q_lim))
            proc.daemon = True
            procs.append(proc)
            proc.start()
        
        [proc.join() for proc in procs]
    
    if mode == 'pool':
        run_pool(interval,
                 path_to_fastq,
                 output,
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
    