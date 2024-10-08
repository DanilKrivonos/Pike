# **Pike: tool for Oxford Nanopore amplicon metagenomics**
<img src="https://github.com/user-attachments/assets/9961d756-2a37-4882-ac1a-b3c1c4df143c" width="250" height="250" align="left">  <br /><br />
Pike is a tool for analyzing amplicon nonopore sequencing data. Pike is capable of analyzing reads of different amplicrone sizes. The algorithm is based on sequential clustering of reads followed by reaching consensus amplicon variants. 
The analysis is possible in two different versions. Single mode allows you to process samples independently, performing clustering and OTU assembly separately for each sample. Pool mode performs joint clustering for all samples, followed by data parsing for individual samples.  
  <br />
  <br />
  
## **Dependencies:**
- numpy>=1.22.4
- scikit-bio>=0.5.9
- biopython>=1.83
- pandas>=2.2.0
- scipy=1.10.1
- scikit-learn>=1.4.1
- hdbscan>=0.8.33
- umap-learn>=0.5.5
- medaka=1.11.3
- mafft
- cutadapt=4.6
- filtlong=0.2.1

## **Installation:**
```
conda create -n pike_env -c bioconda -c conda-forge medaka=1.11.3 cutadapt=4.6 filtlong=0.2.1
conda activate pike_env
pip install pike-meta
```

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
  -k K                  k-mer size number (6 by default)

consensus building arguments:
  -consensus_seq_lim CONSENSUS_SEQ_LIM
                        Minimum number of reads supporting consensus sequence
  -letter_Q_lim LETTER_Q_LIM
                        The minimum quality of the consensus letter. If the quality is less, then N is set (15 by default)
```

### **Single mode:**
```
pike -mode single -fastq example_data/ --trim_primer  -primerF CCTACGGGNGGCWGCAG -primerR GACTACHVGGGTATCTAATCC
```
### **Pool mode:**
```
pike -mode pool -fastq example_data/ --trim_primer  -primerF CCTACGGGNGGCWGCAG -primerR GACTACHVGGGTATCTAATCC
```
### **Taxonomy calling via blastn:**
```
get_taxonomy  -otutab Pike_results/merged_otu_table.tsv -output Pike_results/TAXONOMY -db /path/to/db.fasta
```
### Output directory architecture

```
.
├── merged_otu_table.tsv            summary OTU table for all samples
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
