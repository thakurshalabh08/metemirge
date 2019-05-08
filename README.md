# metemirge
Bioinformatics pipeline for reconstructing marker genes from metagenome sequence data, using EMIRGE.

## Overview

This pipeline is written in snakemake and designed to automate and control the submission of processes to the Synergy server at the University of Calgary in the lab of Dr. Laura Sycuro. Developed by Shalabh Thakur, Laura Sycuro, and Alana Schick.

EMIRGE is a tool for reconstructing full-length marker gene sequences from microbial community short-read sequencing data (Miller et al, 2011, Genome Biology). 

Input: filtered and cleaned fastq files. 

Output: ???

## Installation

To use this pipeline, clone this repository into your project directory using the following command:

```
git clone https://github.com/SycuroLab/metemirge.git
```

Note: you need to have snakemake installed in order to dun this. To install snakemake using conda, run the following line:

```
conda install -c bioconda snakemake
```

See the [snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Home) for details.

## Config file

This pipeline requires a config file, written in yaml, to run. Enter any custom parameters in `config.yaml`. This is the only file that should be modified before running the pipeline. 

The string parameter value should be written within double quotes " ".

The numerical parameter values should be written without an quotes.

| Parameters | Description |
| ---------- | ----------- |
| **General Parameters** |
| list_files | The file path containing list of sample names or ids (string) |
| input_dir | The directory path for input fastq files (string) |
| project_dir | The directory path for the project directory (string) |
| **Cluster Parameters** |
| num_threads | Specify number of threads to use for the analysis. (integer) |
| **Emirge Parameters** |
| fasta_db | Specify full path to the reference fasta file for bowtie (string) |
| bowtie_db | Specify full path to the reference bowtie index file (string) |
| max_read_length | Specify maximum read length for the samples (integer) |
| insert_mean | Specify mean value for the insert size (integer) |
| insert_stddev | Specify standard deviation value for the insert size (decimal) |
| num_iter | Specify number of iteration for EMIRGE as an integer (integer) |
| num_iter_str | Specify number of iteration for EMIRGE as a string within " " (string) |
| **Blast Parameters** |
| blast_db_dir | The base directory path for the blast databases (string) |
| databases | Specify list of blast databases as ["db1", "db2"] (list) |
| evalue | Specify blast evalue cutoff (decimal / scientific) |
| max_target_seqs | Specify maximum target sequences to report in the blast result (integer) |
| outfmt | Specify blast output format. Analysis requires tabular output format. Default: 6 (integer)|
| out_columns | Specify columns to be reported in the blast output. Should not be changed (string) |
| **Blast Parsing Parameters** |
| silva_taxdb | The sqlite database file storing taxonomy information from silva database (string) |
| ncbi_taxdb | The sqlite database file storing taxonomy information from ncbi database (string) |
| percent_identity_threshold | Specify identity threshold for filtering hits (decimal) |
| query_coverage_threhshold | Specify query coverage threshold for filtering hits (decimal) |
| **Select Best-Hit Parameters** |
| db_priority | Specify list of blast databases in order of priority as ["db1", "db2"] (list) |
| num_best_hit | Specify number of best hits to report from each database (integer) |




## Running the pipeline on Synergy

Test the pipeline by running `snakemake -np`. This command prints out the jobs that will be submitted to the Synergy compute cluster without actually submitting them.

To run the pipeline on the cluster, enter the following command from the project directory:

```
snakemake --cluster-config cluster.json --cluster 'bsub -n {cluster.n} -R {cluster.resources} -W {cluster.walllim} -We {cluster.time} -M {cluster.maxmem} -oo {cluster.output} -e {cluster.error}' --jobs 500 --use-conda
```

Note: the file `cluster.json` contains the parameters for the LSF job submission system. They are by default the same for each process, but can be modified in this file.

## Results and log files

All output files for each sample will be placed in the `output/<project-dir>/<sample-dir>` directory and logs of each rule/process will be written to the `logs` directory.

## Pipeline summary

1) Reconstruct marker sequences using EMIRGE.

2) Rename marker sequence file.

3) Blast marker sequence against reference sequence databases.

4) Parse Blast results.

5) Select best blast hit.

6) Prepare Input for Plot.
