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
git clone https://github.com/alanaschick/metemirge.git projectname
```

Note: you need to have snakemake installed in order to dun this. To install snakemake using conda, run the following line:

```
conda install -c bioconda snakemake
```

See the [snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Home) for details.

## Config file

This pipeline requires a config file, written in yaml, to run. Enter any custom parameters in `config.yaml`. This is the only file that should be modified before running the pipeline. 

| Parameters | Description |
| ---------- | ----------- |
| General Parameters |
| list_files | The file containing list of sample names or ids |
| input_dir | The directory path for input fastq files |
| project_dir | The directory path for the project directory |
| Emirge Parameters |

## Running the pipeline on Synergy

Test the pipeline by running `snakemake -np`. This command prints out the jobs that will be submitted to the Synergy compute cluster without actually submitting them.

To run the pipeline on the cluster, enter the following command from the project directory:

```
snakemake --cluster-config cluster.json --cluster 'bsub -n {cluster.n} -R {cluster.resources} -W {cluster.walllim} -We {cluster.time} -M {cluster.maxmem} -oo {cluster.output} -e {cluster.error}' --jobs 500 --use-conda
```

Note: the file `cluster.json` contains the parameters for the LSF job submission system. They are by default the same for each process, but can be modified in this file.

## Results and log files

All output files will be placed in the `results` directory and logs of each rule/process will be written to the `logs` directory.

## Pipeline summary

1) Reconstruct marker sequences using EMIRGE.

2) Rename marker sequence file.

3) Blast marker sequence against reference sequence databases.

4) Parse Blast results.

5) Select best blast hit.

6) Prepare Input for Plot.
