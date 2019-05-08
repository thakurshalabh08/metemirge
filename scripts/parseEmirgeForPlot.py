#!/usr/bin/env python

"""
Shalabh Thakur
May, 2018
Sycuro Lab, University of Calgary
"""
"""
Description: This script reads csv result file from Selected best hits for priority database and parse to create input file with reduced information for generating
taxonomic abundance plot in R
"""
import sys
import os
import glob
import csv
import argparse
import re
import shutil
import logging
import datetime
from argparse import RawTextHelpFormatter
from collections import defaultdict


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
	parser.add_argument("--prefix", help="Name of the sample", required=True)
	parser.add_argument("--input_file", help="Name and path of the input file", required=True)
	parser.add_argument("--output_file", help="name and path of the output file", required=True)

	args=parser.parse_args()

	return vars(args)


def main():

	"""
	Main Function
	"""
	now = datetime.datetime.now()
	
	c_args = parse_args(__file__)

	sample_name=c_args['prefix']
	input_file=os.path.abspath(c_args['input_file'])
	output_file=os.path.abspath(c_args['output_file'])

	#log_file=os.path.join(project_path,sample_name+"_"+project_name.lower()+"_run.log")
	logging.basicConfig(format="[%(asctime)s] %(levelname)-8s %(message)s",datefmt='%a, %d %b %Y %H:%M:%S',level=logging.INFO)

	logging.info("Writing parsed data for plotting results for sample {} in {}\n".format(sample_name,output_file))

	write_plot_file(sample_name,input_file,output_file)



##### FUNCTIONS #######
#### write parsed file for making plots #####
def write_plot_file(sample_name,input_file,output_file):

	
	##### open output file to write ####
	
	plot_data=open(output_file,"w")
	
	### Write Header ###
	plot_data.write("SAMPLE_NAME\tQUERY_ID\tSUBJECT_ID\tSUBJECT_SPECIES\tTAXONOMY\tEMIRGE_NORM_PRIOR\tPERCENT_ABUNDANCE\tIS_SPECIES\n")

	##### read csv file for best-selected blast hits for sample ####

	with open(input_file,"r") as best_blast_hit:

		blast_parsed_dict=csv.DictReader(best_blast_hit,delimiter='	')
		
		for row in blast_parsed_dict:

			query_id=row['QUERY_ID']
			subject_id=row['SUBJECT_ID']
			subject_species=row['SUBJECT_SPECIES']
			taxonomy=row['TAXONOMY']
			emirge_prior=row['EMIRGE_NORM_PRIOR']
			is_species=row['IS_SPECIES']

			percent_abundance=round(float(emirge_prior)*100,2)

			plot_data.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample_name,query_id,subject_id,subject_species,taxonomy,emirge_prior,percent_abundance,is_species))

	plot_data.close()


if __name__ == '__main__':
	main()		
