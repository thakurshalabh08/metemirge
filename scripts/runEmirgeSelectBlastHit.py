#!/usr/bin/env python

"""
Shalabh Thakur
May, 2018
Sycuro Lab, University of Calgary
"""
"""
Description: This script is used to read formatted and parsed BLAST results for EMIRGE sequences and select best hit in order of database priority.
This version of runEmirgeSelectBlastHit.py choose high scoring species level hit from reference databases in order of priority.
If no species level hit found within any priority reference database then choose highest scoring hit from any database.
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
	parser.add_argument("--prefix", default="blast_result", help="prefix to add to the name of output file",required=True)
	parser.add_argument("--input_files", help="Name and path of the parsed tabular blast result file",nargs='+', required=True)
	parser.add_argument("--output_file", help="Name path of the output file", required=True)
	parser.add_argument("--priority_db", help="list of blast database names in order of priority", nargs='+',required=True)
	parser.add_argument("--num_hits", help="Specify number of species-level hits to select for each query from priority datbase", required=True)	
	args=parser.parse_args()

	return vars(args)



##### MAIN FUNCTION #######
def main():

	"""
	Main Function
	"""
	now = datetime.datetime.now()
	
	c_args = parse_args(__file__)

	prefix=c_args['prefix']
	input_files=c_args['input_files']
	output_file=os.path.abspath(c_args['output_file'])
	priority_db_list=c_args['priority_db']
	num_hits=int(c_args['num_hits'])

	logging.basicConfig(format="[%(asctime)s] %(levelname)-8s %(message)s",datefmt='%a, %d %b %Y %H:%M:%S',level=logging.INFO)

	######### Define Reference Blast Database Priority #######
	priority_db_dict={v:v for i,v in enumerate(priority_db_list)}

	#### Read hits from parsed result ########

	db_hits, headers=parse_result(input_files,priority_db_dict)

	write_hits_by_db_priority(db_hits,headers,priority_db_list,num_hits,output_file)

	##### Write parsed result with reduced information for making plots ####

	#os.system("python parseEmirgeForPlot.py --sample_name {} --project_name --project_path {} --input_file {}".format(sample_name,project_name,project_path,output_file))

	

##### FUNCTIONS #######



#### Function to read sample dir names #####

def get_dir_list(directory):
	"""
	Input: directory name
	"""
	return [name for name in os.listdir(directory)
		if os.path.isdir(os.path.join(directory, name))]

#### Function to read sample file names ####
def get_file_name_from_dir(directory):
	"""
	Input: directory name
	"""
	return [name for name in os.listdir(directory)
		if os.path.isfile(os.path.join(directory, name))]



##### Function to parse hits #####

def parse_result(input_files,priority_db_dict):

	"""
	This function read hit from the parsed blast result file into a dictionary
	"""
	
	####Template to search for Parsed Blast Result Dir ####
	db_hits=defaultdict(dict)
	
	headers=None

	##### Loop through each blast result file #########
	for parsed_result_file in input_files:

		parsed_result_file=os.path.abspath(parsed_result_file)

		logging.info("[SELECT_HIT]: Reading parsed blast result file from {}".format(parsed_result_file))

		with open(parsed_result_file,"r") as blast_parsed_tsv:

			blast_parsed_dict=csv.DictReader(blast_parsed_tsv,delimiter='	')

			for line in blast_parsed_dict:

				query_id=line['QUERY_ID']

				db_name=line['BLAST_DB']

				db_name=os.path.basename(db_name)
				db_name=str(os.path.splitext(db_name)[0])

				##### Check if incorrect database name given on command line #######
				if not priority_db_dict[db_name]:

					logging.error("Invalid blast database name. {} database not found in blast result".format(db_name))
					sys.exit()
				else:
					line['BLAST_DB']=db_name

					if query_id in db_hits and db_name in db_hits[query_id]:

						hits=db_hits[query_id][db_name]
						hits.append(line)
						db_hits[query_id][db_name]=hits

					else:
						hits=list()
						hits.append(line)
						db_hits[query_id][db_name]=hits


					if not headers:
						headers = blast_parsed_dict.fieldnames
					
						

	return(db_hits,headers)


###### Function to select hits by database priority #######
def write_hits_by_db_priority(db_hits,headers,priority_db_list,num_hits,output_file):

	"""
	Choose top N hits at species level based on database priority. If hits do not belong to species level
	then choose hit based on highest bitscore from any database
	"""

	logging.info("Saving selected top hits in: {}".format(output_file))

	output_tsv_handle = open(output_file,"w")
	tsv_writer = csv.DictWriter(output_tsv_handle, headers, dialect='excel-tab')
	tsv_writer.writeheader()


	for query_id in db_hits:

		##### Sort match for each query id by highest bitscore #####

		query_db_dict=db_hits[query_id]

		priority_hits=list()
		no_priority_hits=list()

		for db_name in priority_db_list:

			query_line_list=query_db_dict[db_name]
		
			for query_line in query_line_list:

				is_species=str(query_line['IS_SPECIES'])
			
				##### if query line has species level hit in priority database, 
				##### use highest bitscore hit only
				if is_species=="True":
					priority_hits.append(query_line)
					break
				else:
				#### if there is no species level-hit in priority db or
				#### species-level hit in non-priority database, use highest bit score hit
					no_priority_hits.append(query_line)

			#### If species-level hit found in priority order, print and exit the loop ####		
			if len(priority_hits)>=1:	
				tsv_writer.writerow(priority_hits[0])
				break

		##### print [num_hits] hits with highest bit score from any database when there is no species-level hit from any priority database ####
		if len(priority_hits)==0:

			highscore_hit=sorted(no_priority_hits, key=lambda k: k['TOTAL_BITSCORE'], reverse=True)

			for index, hit in enumerate(highscore_hit):
				if index < num_hits: 
					#print "{}\n\n".format(hit)
					tsv_writer.writerow(hit)
				

if __name__ == '__main__':
	main()



