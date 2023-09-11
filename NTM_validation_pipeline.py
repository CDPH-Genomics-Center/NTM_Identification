#!/usr/bin/env python3

"""
Description
"""

#Imported modules
import os
import time
import json
import shutil
import argparse
from math import ceil
from statistics import mean
from subprocess import call, run
from multiprocessing import cpu_count


#Authorship information
__author__ = "Angus Angermeyer"
__copyright__ = "Copyright 2023, Dr. Angus Angermeyer"
__license__ = "GPL"
__version__ = "0.1.0"
__email__ = "angus.angermeyer@gmail.com"



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
thread_count = 4




def log_result(item, value, font=None):
	offset = 7
	gap = '\t' * (offset-ceil(len(item)/4+0.01))
	
# 	
# 	if len(item) < 4:
# 		gap = "\t\t\t\t\t\t"
# 	
# 	elif len(item) > 3:
# 		gap = "\t\t\t\t\t"
# 
# 	elif len(item) > 7:
# 		gap = "\t\t\t\t"
# 	
# 	elif len(item) > 11:
# 		gap = "\t\t\t"
# 	
# 	elif len(item) > 15:
# 		gap = "\t\t"
# 	
# 	elif len(item)  > 19:
# 		gap = "\t"
# 	
# 	else:
# 		gap = "\t\t\t\t\t"
	
	with open(result_log_name, 'a') as file:
		if font == 'header':
# 			file.write('#'*50)
# 			file.write(f"\n{item}\n\n")
			file.write(f"\n{'#'*50}\n{item}\n\n")
		else:	
			file.write(str(item) + gap + str(value) +'\n')


def fastp_filtering(input_reads, min_fastq_length, min_fastq_qscore):
	"""
	Description
	"""
	
	run(f"fastp -i {input_reads} -o trimmed.fastq.gz -l {min_fastq_length} -q {min_fastq_qscore} -w {args.threads} >> {console_log_name} 2>&1", shell=True)


	with open('fastp.json', 'r') as fastp_file:
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


def flye_assembly(long_reads):
	"""
	Description
	"""
	
	run("flye "
		f"--nano-corr {long_reads} "
		f"--out-dir flye_assembly "
		f"--threads {thread_count} "
		f">> {console_log_name} 2>&1",
		shell=True
	)
	
	with open("flye_assembly/assembly_info.txt", 'r') as assembly_info:
		assembly_info.readline()
		
		contig_lens = []
		contig_covs = []
		
		for line in assembly_info:
			line = line.split()
			contig_lens.append(int(line[1]))
			contig_covs.append(int(line[2]))
			
		mean_len = int(mean(contig_lens))
		mean_cov = round(mean(contig_covs),1)
		total_len = sum(contig_lens)
		N50 = calculate_n50(contig_lens)
		
		log_result("Flye version", "2.9.2")
		log_result("Mean contig length", mean_len)
		log_result("Mean contig coverage", mean_cov)
		log_result("Assembly length", total_len)
		log_result("N50", N50)
		
				
	os.rename("flye_assembly/assembly.fasta", f"{run_name}_assembly.fasta")
	shutil.rmtree("flye_assembly")

def calculate_n50(contig_lengths):
	contig_lengths.sort(reverse=True)
	total_length = sum(contig_lengths)
	half_total_length = total_length / 2

	return next(length for length in contig_lengths if (half_total_length := half_total_length - length) <= 0)


# def unicycler_assembly(long_reads, min_fasta_length):
# 	"""
# 	Description
# 	"""
# 
# 	run("unicycler "
# 		f"--out unicycler_result " 
# 		f"--long {long_reads} "
# 		f"--min_fasta_length {min_fasta_length} "
# 		"--verbosity 2 "
# 		f">> {console_log_name} 2>&1",
# 		shell=True
# 		)
# 		
# 		
# 	os.rename("unicycler_result/assembly.fasta", f"{run_name}_assembly.fasta")
# 	shutil.rmtree("unicycler_result")	


def fastani_compute(assembly, reference_folder="ANI_references"):
	"""
	Description
	"""
	
	with open("reference_list.txt", 'w') as reference_list_file:
		ref_files = os.listdir(reference_folder)
		for ref_file in ref_files:
			reference_list_file.write(f"{reference_folder}/{ref_file}\n")

	
	run("fastANI "
		f"-q {assembly} "
		f"--rl reference_list.txt "
		f"-o ANI_result.txt "
		"--threads 4 "
		f">> {console_log_name} 2>&1",
		shell=True)
	
	
	with open("ANI_result.txt", 'r') as ANI_result:
		top_hit = ANI_result.readline().split()
		top_ref = top_hit[1]
		top_ANI = top_hit[2]
		
		with open(top_ref, 'r') as ref_assembly:
			header = ref_assembly.readline()
			accession = header[1:header.index(' ')]
			species = header[header.index(' ')+1:header.index(' strain')]

		log_result("fastANI version", "1.33")
		log_result("Species", species)
		log_result("ANI", top_ANI)
		log_result("Accession", accession)
	
	os.remove("reference_list.txt")
	os.remove("ANI_result.txt")





log_result("NTM Validation Pipeline (pre-alpha)", None, "header")
log_result("Sample name", args.sample_name)
log_result("Date", run_name.split('_')[-2])
log_result("Time", run_name.split('_')[-1].replace('-',':'))




log_result("Read QC and filtering with fastp", None, "header")
trimmed_reads = fastp_filtering(args.nanopore_reads, args.min_fastq_length, args.min_fastq_qscore)

log_result("De novo assembley with Flye", None, "header")
flye_assembly(trimmed_reads)

log_result("Species ANI identification with fastANI", None, "header")
fastani_compute(f"{run_name}_assembly.fasta")
























