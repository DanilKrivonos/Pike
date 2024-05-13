import os
import numpy as np
from os import listdir
from umap import UMAP
from umap.parametric_umap import ParametricUMAP
from hdbscan import HDBSCAN
from subprocess import call
from skbio.stats.composition import clr
from pandas import DataFrame
from multiprocessing import Process
from sklearn.decomposition import PCA
from src.get_pool import create_pool_dirs, collect_pool_features, run_pool
from src.get_features import collect_features
from src.get_visualization import get_visualisation
from src.get_filtered_clusters import filter_cluster
from src.get_good_reads import run_trimming, run_filtering
from src.get_consensus import get_msa_info, get_consensus, medaka_run, prepare_output
import tensorflow as tf

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
        barcode_path = path_to_fastq

        #try:
        create_dirs(output, barcode) #creating and removing directoryes
        
        print(f' R E A D S   P R E P R O C E S S I N G   F O R    :   {barcode} ')
        print('===================================================================')

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

        print(f' F E A T U R E S   C O L L E C T I N G   F O R    :   {barcode} ')
        print('===================================================================')

        K_MERS_FREQ, GC_CONTENT, READ_ID, READ_Seq, READ_Q, LENS, QUALITY, BARCODE_ID = collect_features(barcode_path, 
                                                                                            barcode, 
                                                                                            usereads,
                                                                                            k, 
                                                                                            read_q_score) #Features collection
                    
        #___Staged decomposition_________________________________________________________________________________________________________
        print(f' C L U S T E R I N G   S T A G E   F O R          :   {barcode} ')
        print('===================================================================')
        K_MERS_FREQ = np.array(K_MERS_FREQ)
        idx = np.argwhere(np.all(K_MERS_FREQ[..., :] == 0, axis=0))
        K_MERS_FREQ = np.delete(K_MERS_FREQ, idx, axis=1)
        K_MERS_FREQ += 1
        normalize = lambda x: x / np.sum(x)
        K_MERS_FREQ = np.array(list(map(normalize, K_MERS_FREQ)))
        clr_data = clr(K_MERS_FREQ)
        import scipy
        #clr_data = clr_data
        pca_model = PCA(n_components=15,
                        random_state=0)
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
                        'QUALITY' : QUALITY,
                        'K-mers signature' : list(K_MERS_FREQ)}
        RESULT_DF = DataFrame(RESULT_DICT) #All information mearged in df

        #___Cluster identification________________________________________________________________________________________________________
        hdbscan = HDBSCAN(min_cluster_size=hdbscan_neighbours,
                        cluster_selection_epsilon=0.1, 
                        gen_min_span_tree=True,
                        metric='euclidean')
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

                read_record = f'@{filtered_data["Read ID"][idx]}\n{filtered_data["READ_Seq"][idx]}\n+\n{filtered_data["READ_Q"][idx]}\n'
                fasta_record = f'>{filtered_data["Read ID"][idx]}\n{filtered_data["READ_Seq"][idx]}\n'
                clt_fastq.write(read_record) #fastq recording
                quality_dict[filtered_data["Read ID"][idx]] = [i for i in RESULT_DF['READ_Q'][idx]] #letter quality adding
                clt_fasta.write(fasta_record) #fasta recording
            
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
        print(f' P O L I S H   C O N S E N S U S    F O R         :   {barcode} ')
        print('===================================================================')

        medaka_run(output, barcode) # run medaka polishing

        print(f' O U T P U T   P R E P A R A T I O N    F O R     :   {barcode} ')
        print('===================================================================')

        prepare_output(output, barcode)

      #  except:
      #      print(f'ERROR WITH {barcode}')

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
        
        create_pool_dirs(output)

        for interval in intervals:
            
            proc = Process(target=collect_pool_features, 
                        args=(interval,
                              path_to_fastq,
                              output,
                              usereads,
                              k,
                              trim_primer, 
                              primerF, 
                              primerR,
                              minlen, 
                              maxlen,
                              read_q_score))
            proc.daemon = True
            procs.append(proc)
            proc.start()
        
        [proc.join() for proc in procs]
            
        run_pool(output,
                 threads,
                 umap_neighbours,
                 hdbscan_neighbours,
                 visualize,
                 consensus_seq_lim,
                 letter_Q_lim)
        
