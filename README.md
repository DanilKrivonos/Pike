# Pike: metagenomic tool for noisy reads analyzing
Pike is a tool for analyzing noisy reads. The algorithm is based on the clustering of reads with the subsequent construction of their multiple pulling out, as a result of which a consensus sequence is built.


## **Dependencies:**
- pandas
- biopython
- numpy
- sklearn
- umap-learn
- hdbscan
- medaka
- mafft

## **Installation:**

## **Ussages:**
### Parameters
```
P I K E : metagenomic tool for noisy reads analyzing

optional arguments:
  -h, --help            show this help message and exit

main arguments:
  -fastq FASTQ          Directory with merged read files
  -output OUTPUT        Directory with output files

advance arguments:
  -mode MODE            Choice of analysis mode. single - each sample is analyzed separately, pool - joint analysis of all samples (single by default)
  -usereads USEREADS    Maximum number of reads per sample for analysis (Everything is used by default)
  -threads THREADS      Number of threads

read quality arguments:
  --trim_primer         Trim primers (If this mod is selected, then it is necessary to apply primer sequences)
  -primerF PRIMERF      Forward primer sequence
  -primerR PRIMERR      Reverse primer sequence
  -read_q_score READ_Q_SCORE
                        Median value of read quality
  -minlen MINLEN        Minimum expected amplicon length (350 by default)
  -maxlen MAXLEN        Maximum expected amplicon length (600 by default)

clustering arguments:
  -umap_neighbours UMAP_NEIGHBOURS
                        Number of nearest neighbors for UMAP (30 by default)
  -cluster_size CLUSTER_SIZE
                        Number of dots for HDBSCAN (30 by default)
  --visualize           Create clustering figure (false by default)

consensus building arguments:
  -consensus_seq_lim CONSENSUS_SEQ_LIM
                        Minimum number of reads supporting consensus sequence
  -letter_Q_lim LETTER_Q_LIM
                        The minimum quality of the consensus letter. If the quality is less, then N is set (15 by default)
```

### **Single mode:**

### **Pool mode:**



### Output directory architecture

```
.
├── read_preprocessing              directory with read preprocessing results
|   ├── cutadapt_round1             directory with cutadapt round 1 results (cut F primer)
|   ├── cutadapt_round2             directory with cutadapt round 2 results (cut R primer)
│   └── filtered_reads              directory with filtered reads
|
├── work_dir                        directory with working files
│   ├── barcode01
│   |   ├── clusters_data_fastq     reads divided into clusters
│   |   ├── clusters_data_fasta     fasta for reads broken into clusters
│   |   ├── msa                     multiple alignment for reads broken into clusters
│   |   ├── proto_consensus         consensus built from multiple alignment
│   |   ├── medaka                  Medaka consensus results
│   |   └── metadata_df.tsv         table with meta information on reads
│  ...
|   └── barcodeN 
└── results                         *Pike results*
    ├── barcode01
    |   ├── results.tsv             table with collected consensus and their representations
    |   └── barcode01.pdf           clustering visualization (by request --visualize)
   ...
    └── barcodeN 
```

## Reference
Add later 