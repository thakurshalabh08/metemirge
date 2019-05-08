# ************************************
# * Snakefile for metemirge pipeline *
# ************************************

# **** Variables ****

configfile: "config.yaml"

# Specify the list of files to run the pipeline on
import os
import re
import pandas as pd

SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Rules ****

rule all:
    input: expand(config["project_dir"]+"emirge/{sample}/{sample}_emirge.fasta",sample=SAMPLES),
           expand(config["project_dir"]+"emirge/{sample}/blast/{database}/{sample}_raw_blast_table.tsv",sample=SAMPLES,database=config["database"]),
           expand(config["project_dir"]+"emirge/{sample}/blast/{database}/{sample}_parsed_blast_table.tsv",sample=SAMPLES,database=config["database"]),
           expand(config["project_dir"]+"emirge/{sample}/{sample}_best_hit.tsv",sample=SAMPLES),
	       expand(config["project_dir"]+"emirge/{sample}/{sample}_abundance_plot_data.tsv",sample=SAMPLES),
           config["project_dir"]+"emirge_plot/abundance_plot.pdf"

rule runemirge:
    input:
        r1 = config["input_dir"]+"{sample}_read1.fastq",
        r2 = config["input_dir"]+"{sample}_read2.fastq"
    output:
        out1 = config["project_dir"]+"emirge/{sample}/iter."+config["num_iter_str"]+"priors.initialized.txt",
        out2 = directory(config["project_dir"]+"emirge/{sample}/iter."+config["num_iter_str"])
    conda: "envs/emirge_env.yaml"
    params:
        outdir = config["project_dir"]+"emirge/{sample}/"
    shell:
        "emirge.py {params.outdir} -1 {input.r1} -2 {input.r2} "
        "-f {config[fasta_db]} -b {config[bowtie_db]} -l {config[max_read_length]} "
        "-i {config[insert_mean]} -s {config[insert_stddev]} -n {config[num_iter]} "
        "-a {config[num_threads]} --phred33"

rule rename:
    input: config["project_dir"]+"emirge/{sample}/iter."+config["num_iter_str"]
    output: config["project_dir"]+"emirge/{sample}/{sample}_emirge.fasta"
    conda: "envs/emirge_env.yaml"
    shell: "emirge_rename_fasta.py {input} > {output}"

rule blast:
    input: config["project_dir"]+"emirge/{sample}/{sample}_emirge.fasta"
    output: expand(config["project_dir"]+"emirge/{{sample}}/blast/{database}/{{sample}}_raw_blast_table.tsv",database = config["database"])
    conda: "envs/blast_env.yaml"
    params:
         out_dir=config["project_dir"]+"emirge/{sample}/blast"
    shell: 
         "python scripts/runBlast.py --query_file {input} --db_dir {config[blast_db_dir]} "
         "--db_names {config[database]} --output_dir {params.out_dir}  --evalue {config[evalue]} "
         "--max_target_seq {config[max_target_seqs]} --output_format {config[outfmt]} "
         "--output_columns {config[out_columns]} --prefix {wildcards.sample}"

rule parseblast:
    input: 
          blast_raw = expand(config["project_dir"]+"emirge/{{sample}}/blast/{database}/{{sample}}_raw_blast_table.tsv",database = config["database"]),
          emirge_fasta = config["project_dir"]+"emirge/{sample}/{sample}_emirge.fasta"
    output: expand(config["project_dir"]+"emirge/{{sample}}/blast/{database}/{{sample}}_parsed_blast_table.tsv",database = config["database"])
    conda: "envs/biopy_env.yaml"
    shell: "python scripts/runEmirgeParseBlast.py --prefix {wildcards.sample} --blast_raw {input.blast_raw} "
           "--emirge_file {input.emirge_fasta} --identity_cutoff {config[percent_identity_threshold]} "
           "--qcov_cutoff {config[query_coverage_threshold]} --silva_taxadb {config[silva_taxdb]} "
           "--ncbi_taxadb {config[ncbi_taxdb]}"

rule selectbest:
    input: expand(config["project_dir"]+"emirge/{{sample}}/blast/{database}/{{sample}}_parsed_blast_table.tsv",database = config["database"])
    output: config["project_dir"]+"emirge/{sample}/{sample}_best_hit.tsv"
    conda: "envs/biopy_env.yaml"
    shell: "python scripts/runEmirgeSelectBlastHit.py --prefix {wildcards.sample} --input_files {input} --output_file {output} "
           "--priority_db {config[db_priority]} --num_hits {config[num_best_hit]}"

rule inputforplot:
    input: config["project_dir"]+"emirge/{sample}/{sample}_best_hit.tsv"
    output: config["project_dir"]+"emirge/{sample}/{sample}_abundance_plot_data.tsv"
    conda: "envs/biopy_env.yaml"
    shell: "python scripts/parseEmirgeForPlot.py --prefix {wildcards.sample} --input_file {input} --output_file {output}"

#rule plotabundance:
#    input: directory(config["project_dir"]+"emirge")
#    output: config["project_dir"]+"emirge_plot/abundance_plot.pdf"
#    conda: "envs/r_env.yaml"
#    shell: "Rscript scripts/plotEmirgeAbundance.R {input} {output}"

