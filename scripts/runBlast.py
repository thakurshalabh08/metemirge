#!/usr/bin/env python
"""
Author: Shalabh Thakur
Description: Wrapper to run blast through multiple databases
"""
import os
import sys
import argparse
import logging
import pandas as pd
import numpy as np
from argparse import RawTextHelpFormatter

def parse_args(desc):

	"""
	Command line option parser
	
	:param desc: Short program description.
	:type desc: str
	
	:return arg: dict of command line arguments with key=command line
        argument and val=argument value
        
	:return_type: dict
	"""

	parser = argparse.ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter)
	parser._optionals.title="Arguments"
	parser.add_argument("--query_file",help="name of the sample", required=True)
	parser.add_argument("--db_names", help="list of database name for blast", nargs='+',required=True)
	parser.add_argument("--db_dir", help="root path for all blast database", required=True)
	parser.add_argument("--evalue", type=float, help="E-value threshold for BLAST", required=True)
	parser.add_argument("--max_target_seqs", type=int, help="Maximum number of target sequences to be reported by BLAST", required=True)
	parser.add_argument("--output_dir", help="path to blast output directory", required=True)
	parser.add_argument("--output_format", type=int, help="Specify format from BLAST result", required=True)
	parser.add_argument("--output_columns", help="specify columns to be reported in BLAST output", nargs='+',required=True)
	parser.add_argument("--prefix", default="blast_result", help="prefix to add to the name of blast output file")

	
	args=parser.parse_args()

	return vars(args)

##### MAIN FUNCTION #######
def main():

	"""
	Main Function
	"""
	
	c_args = parse_args(__file__)

	logging.basicConfig(format="[%(asctime)s] %(levelname)-8s %(message)s",datefmt='%a, %d %b %Y %H:%M:%S',level=logging.INFO)

	query_file=os.path.abspath(c_args["query_file"])
	db_dir=os.path.abspath(c_args["db_dir"])
	out_dir=os.path.abspath(c_args["output_dir"])
	db_list=c_args["db_names"]
	evalue=float(c_args["evalue"])
	max_target=int(c_args["max_target_seqs"])
	out_format=c_args["output_format"]
	out_columns=" ".join(c_args["output_columns"])
	prefix=c_args["prefix"]


	#### run BLAST through each database in the list ###

	logging.info("Running BLAST:")

	for db_name in db_list:

		if not os.path.exists(os.path.join(out_dir,db_name)):
			os.makedirs(os.path.join(out_dir,db_name))
			logging.info("Creating output folder for database = "+db_name)
			
		db_file=os.path.join(db_dir,db_name,db_name)
		out_file=os.path.join(out_dir,db_name,prefix+"_raw_blast_table.tsv")

		logging.info("query = "+query_file)
		logging.info("db = "+db_file)
		logging.info("out = "+out_file)

		os.system("blastn -query {query_file} -db {db_file} -out {out_file} -evalue {evalue} -max_target_seqs {max_target} -outfmt '{out_format} {out_columns}'".
			format(query_file=query_file,db_file=db_file,out_file=out_file,evalue=evalue,max_target=max_target,out_format=out_format,out_columns=out_columns))

		### adding column to blast result with database name ###
		header_list=c_args["output_columns"]

		df=pd.read_csv(out_file, sep="\t",header=None)
		df.columns=header_list
		df["database"]=db_name
		df.replace('',np.nan,inplace=True)
		df.to_csv(out_file,sep="\t",index=False)


	sys.exit()


if __name__ == '__main__':
	main()
