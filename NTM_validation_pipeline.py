#!/usr/bin/env python3

"""

"""

#Imported modules
import os
import time
import json
import argparse
from subprocess import call
from multiprocessing import cpu_count


#Authorship information
__author__ = "Angus Angermeyer"
__copyright__ = "Copyright 2023, Dr. Angus Angermeyer"
__license__ = "GPL"
__version__ = "0.1.0"
__email__ = "angus.angermeyer@gmail.com"


f"NTM_{time.strftime('%d-%m-%y_%H-%M')}"

parser = argparse.ArgumentParser(
	description = __doc__, usage="\nPlease see required inputs below. Use '-h' for details.")

parser.add_argument('-s', '--sample_name', type=str, required=True,
	help = "Provide a uniqure sample name for this run.")

parser.add_argument('-n', '--nanopore_reads', type=str, required=True,
	help = "Nanopore read file. Please concatenate multiple reads in to one file. 'cat *.fastq.gz > merged.fq.gz'")

parser.add_argument('--min_fastq_length', type=str, default=500,
	help = "Optional. Default = 500")

parser.add_argument('--min_fastq_qscore', type=str, default=7,
	help = "Optional. Default = 7")

parser.add_argument('--min_contig_length', type=str, default=500,
	help = "Optional. Default = 500")

parser.add_argument('-t', '--threads', type=str, default=cpu_count(),
	help="Optional. Default = Max. available") 

args = parser.parse_args()

run_name = args.sample_name +"_"+ time.strftime('%d-%m-%y_%H-%M')
result_log_name = run_name + "_results_log.txt"
console_log_name = run_name + "_console_log.txt"


def log_result(item, value):
	with open(result_log_name, 'a') as file:
		file.write(str(item) +'\t'+ str(value) +'\n')


def fastp_filtering(input_reads, min_fastq_length, min_fastq_qscore):
	"""
	test
	"""
	call(f"fastp -i {input_reads} -o trimmed.fastq.gz -l {min_fastq_length} -q {min_fastq_qscore} -w {args.threads} >> {console_log_name} 2>&1", shell=True)



	with open('fastp.json') as fastp_file:
		fastp_json = json.load(fastp_file)
		
		log_result("fastp version", fastp_json['summary']['fastp_version'])
		log_result("Total reads (raw)", fastp_json['summary']['before_filtering']['total_reads'])
		log_result("Total bases (raw)", fastp_json['summary']['before_filtering']['total_bases'])
		log_result("Mean length (raw)", fastp_json['summary']['before_filtering']['read1_mean_length'])
		log_result("GC content (raw)", fastp_json['summary']['before_filtering']['gc_content'])
		
		log_result("Total reads (filtered)", fastp_json['summary']['after_filtering']['total_reads'])
		log_result("Total bases (filtered)", fastp_json['summary']['after_filtering']['total_bases'])
		log_result("Mean length (filtered)", fastp_json['summary']['after_filtering']['read1_mean_length'])
		log_result("GC content (filtered)", fastp_json['summary']['after_filtering']['gc_content'])
		
	os.remove("fastp.html")
	os.remove("fastp.json")


	if os.path.isfile("trimmed.fastq.gz"):
		return("trimmed.fastq.gz")
	else:
		log_result("FAILURE", "fastp filtering did not produce trimmed reads file.")







log_result("Project", run_name)



trimmed_reads = fastp_filtering(args.nanopore_reads, args.min_fastq_length, args.min_fastq_qscore)






# ref_list = "paths-to-ref-files.txt"





