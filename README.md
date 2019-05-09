# metemirge
Bioinformatics pipeline for reconstructing marker genes from metagenome sequence data, using EMIRGE.

## Overview

This pipeline is written in snakemake and designed to automate and control the submission of processes to the Synergy server at the University of Calgary in the lab of Dr. Laura Sycuro. Developed by Shalabh Thakur, Laura Sycuro, and Alana Schick.

EMIRGE is a tool for reconstructing full-length marker gene sequences from microbial community short-read sequencing data (Miller et al, 2011, Genome Biology). 

Input: filtered and cleaned fastq files. 

Output: 
1) Reconstructed marker sequences
2) Taxonomy assignment
3) Taxonomic Relative abundance

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
| **Best-Hit Parameters** |
| db_priority | Specify list of blast databases in order of priority as ["db1", "db2"] (list) |
| num_best_hit | Specify number of best hits to report from each database (integer) |


## Scripts folder

 This folder contains customized scripts to run within the pipeline at each step.


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

**1) Reconstruct marker sequences using EMIRGE.**

This step runs EMIRGE program with `emirge.py` script for reconstructing marker sequences using reference sequence database.

```
emirge.py {outdir} -1 {read1} -2 {read2} -f {fasta_db} -b {bowtie_db} -l {max_read_length} -i {insert_mean} -s {insert_stddev} -n {num_iter} -a {num_threads} --phred33
```      

**2) Rename marker sequence file.**

This step renames the reconstructed marker sequence file from final EMIRGE iteration using `emirge_rename_fasta.py`.

```
emirge_rename_fasta.py {path_to_final_iteration_dir} > {renamed_output_file}
```

**3) Blast marker sequence against reference sequence databases.**

This step perform BLAST for reconstructed marker sequences against reference blast databases.

```
python scripts/runBlast.py --prefix {sample_name} --query_file {query_file} --db_dir {blast_db_dir} --db_names {database_list} --output_dir {out_dir}  --evalue {evalue} --max_target_seq {max_target_seqs} --output_format 6 --output_columns "qseqid qacc qlen sseqid sacc slen stitle qstart qend sstart send length evalue bitscore pident qcovs nident mismatch positive gaps qframe sframe staxids sskingdoms"
```

**4) Parse Blast results.**

This step parse tabular BLAST result to filter top-scoring hits within each database.

```
python scripts/runEmirgeParseBlast.py --prefix {sample_name} --blast_raw {blast_raw_result} --emirge_file {emirge_renamed_fasta_file} --identity_cutoff {percent_identity_threshold} --qcov_cutoff {query_coverage_threshold} --silva_taxadb {silva_taxdb} --ncbi_taxadb {ncbi_taxdb}
```

**5) Select best blast hit.**

This step select best BLAST hit for each query from given list of blast databases based on specified database priority.

```
python scripts/runEmirgeSelectBlastHit.py --prefix {sample_name} --input_files {parsed_blast_result_files} --output_file {output_file} --priority_db {db_priority} --num_hits {num_best_hit}
```

**6) Prepare Input for Plot.**

This step parse final best BLAST hit file for each sample and prepare input file for plotting taxonomic abundance in R.

```
python scripts/parseEmirgeForPlot.py --prefix {sample_name} --input_file {input_file} --output_file {output_file}
```
