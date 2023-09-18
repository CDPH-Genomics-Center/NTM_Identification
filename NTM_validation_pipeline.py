#!/usr/bin/env python3

"""
Description
"""

#Imported modules
import os
import glob
import time
import json
import shutil
import argparse
import newick_tree_plotter
from math import ceil
from pathlib import Path
from statistics import mean
from subprocess import run, check_output, STDOUT
from multiprocessing import cpu_count


#Authorship information
__author__ = "Angus Angermeyer"
__copyright__ = "Copyright 2023, Dr. Angus Angermeyer (CDPH)"
__license__ = "GPL"
__version__ = "0.1.0"
__email__ = "angus.angermeyer@gmail.com"



parser = argparse.ArgumentParser(
	description = __doc__, usage="\nPlease see required inputs below. Use '-h' for details.")

parser.add_argument('-s', '--sample_name', type=str, required=True,
	help = "Provide a uniqure sample name for this run.")

parser.add_argument('-n', '--nanopore_reads', type=str, required=True,
	help = "Nanopore read file. Please concatenate multiple reads in to one file. 'cat *.fastq.gz > merged.fq.gz'")

parser.add_argument('--NTM_ref_dir', type=str, default="NTM_references/",
	help = "")

parser.add_argument('--target_seq_dir', type=str, default="target_sequences/",
	help = "")

parser.add_argument('--min_fastq_length', type=str, default=500,
	help = "Optional. Default = 500")

parser.add_argument('--min_fastq_qscore', type=str, default=7,
	help = "Optional. Default = 7")

parser.add_argument('--min_contig_length', type=str, default=500,
	help = "Optional. Default = 500")

parser.add_argument('--blast_coverage', type=str, default=90,
	help = "Optional. Default = 90")

parser.add_argument('--bootstrap_threshold', type=str, default=70,
	help = "Optional. Default = 70")

parser.add_argument('-t', '--threads', type=str, default=cpu_count(),
	help="Optional. Default = Max. available") 



args = parser.parse_args()
run_name = args.sample_name +"_"+ time.strftime('%d-%m-%y_%H-%M')
result_log_name = run_name + "_results_log.txt"
console_log_name = run_name + "_console_log.txt"
thread_count = 4



def log_result(item, value, font=None):
	#This is just some dumb logic to make the results file look pretty.
	#Adds tabs to mae everything line up. Increase offset to make space for longer item names.
	offset = 8
	gap = '\t' * (offset-ceil(len(item.lstrip('\n'))/4+0.01))
	
	#Different type of results logging.
	with open(result_log_name, 'a') as file:
		if font == 'header':
			file.write(f"\n{'#'*50}\n{item}\n\n")

		elif font == 'version':
			terminal_out = check_output(value.split(), text=True, stderr=STDOUT)
			version_line = terminal_out.split('\n')[0]			
			file.write(str(item) + gap + version_line +'\n')

		else:	
			file.write(str(item) + gap + str(value) +'\n')


def fastp_filtering(input_reads, min_fastq_length, min_fastq_qscore):
	"""
	Description
	"""
	
	run(f"fastp -i {input_reads} -o trimmed.fastq.gz -l {min_fastq_length} -q {min_fastq_qscore} -w {args.threads} >> {console_log_name} 2>&1", shell=True)


	with open('fastp.json', 'r') as fastp_file:
		fastp_json = json.load(fastp_file)
		
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


def fastani_compute(assembly, reference_folder=args.NTM_ref_dir):
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

		log_result("Species", species)
		log_result("ANI", top_ANI)
		log_result("Accession", accession)
	
	os.remove("reference_list.txt")
	os.remove("ANI_result.txt")


def extract_target(sample_genome, target_dir, reference_dir, qcov_hsp_perc=90):
	"""
	Description
	"""
	
	from Bio import SeqIO
	from subprocess import check_output
	
	
	def blast_top_hit(subject_file, query_file, coverage_cutoff):
		
		blast_call = ("blastn "
		f"-subject {subject_file} "
		f"-query {query_file} "
		f"-outfmt 6 "
		f"-qcov_hsp_perc {coverage_cutoff} ")
# 		f">> {console_log_name} 2>&1" ) #Some shell parsing issue with this. Maybe should output results to file?

		blast_call = blast_call.split()
		blast_output = check_output(blast_call, text=True)
		
		
	#	query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
		blast_output = blast_output.split("\n")

		hit_list = []
		for entry in blast_output:
			entry = entry.split("\t")		
			hit_list.append([entry[-1],entry]) 
	
		hit_list.sort()	
		top_hit = hit_list[-1][1]
		
		if len(top_hit) > 1:
# 			print(top_hit)
		
			for seq_record in SeqIO.parse(subject_file, "fasta"):
				if seq_record.id == top_hit[1]:
					sstrt = int(top_hit[8])
					sstop = int(top_hit[9])
				
					if sstrt < sstop:
						sequence = str(seq_record.seq[sstrt-1:sstop]) #Weird slicing accounts for indexing differential b/w blast and biopython

					elif sstrt > sstop:
						sequence = str(seq_record.seq[sstop-1:sstrt].reverse_complement())
					
					hit_accession = top_hit[0]
					
					return(hit_accession, sequence)
		
		else:
			return None


	#Extract targets from sample genome#
	
	sample_name = Path(sample_genome).stem
	
	log_result("\nTarget:", "Ref accession:")
	
	for target_seq in Path(target_dir).glob("*.fasta"):
		target_name = target_seq.stem
		
		hit_accession, sequence = blast_top_hit(sample_genome, target_seq, qcov_hsp_perc)
		

		with open(f"{sample_name}_{target_name}.fasta", 'w') as fasta_output:
			fasta_output.write(f">{args.sample_name}_{target_name}\n") ##Grabbing sample name from global here. Should fix.
			fasta_output.write(sequence+"\n")

			log_result(f"{target_name}", f"{hit_accession}")


		#Then ensure that all reference genomes in reference_dir have been blasted for target#
	
		if os.path.isdir(f'{reference_dir}/{target_name}') != True:
			os.mkdir(f'{reference_dir}/{target_name}')
			
			
		for reference_seq in Path(reference_dir).glob("*.fna"):
			ref_name = reference_seq.stem
			
			if os.path.isfile(f'{reference_dir}/{target_name}/{ref_name}_{target_name}.fasta') == True:
				pass
			
			else:
				try:
					ref_hit_accession, ref_sequence = blast_top_hit(reference_seq, target_seq, qcov_hsp_perc)
	
					with open(f"{reference_dir}/{target_name}/{ref_name}_{target_name}.fasta", 'w') as fasta_output:
						ref_record = next(SeqIO.parse(reference_seq, 'fasta'))
						ref_header_list = ref_record.description.split(' ')
						genus = ref_header_list[1]
						species = ref_header_list[2]
						accession = ref_header_list[0]
						
						if genus == 'Mycobacterium':
							genus = 'M.'
						
						fasta_output.write(f">{genus}_{species}_{accession}\n")
						fasta_output.write(ref_sequence+"\n")

				except TypeError: #If no blast hit, it will return nonetype instead of sequence
					with open(console_log_name, 'a') as log:
						log.write(f"\n\nNo blast hit:\t{reference_seq}\t{target_seq}\n\n")
					



def concatenate_align(sample_target_seq, ref_target_dir):
	"""
	Description
	"""
	
	sample_target_name = Path(sample_target_seq).stem

	with open(f"{sample_target_name}_references.fasta", 'w') as combined_fasta:
		with open(sample_target_seq, 'r') as sample_fasta:
			combined_fasta.write(sample_fasta.read())
		
		ref_count = 1
		for ref_target_seq in Path(ref_target_dir).glob("*.fasta"):
			with open(ref_target_seq, 'r') as ref_target_fasta:
				combined_fasta.write(ref_target_fasta.read())
			ref_count+=1

	run("muscle "
		f"-in {sample_target_name}_references.fasta "
		f"-out {sample_target_name}_references.aln "
		"-quiet "
		f">> {console_log_name} 2>&1",
		shell=True)

	log_result(Path(ref_target_dir).name, f"{ref_count} aligned seqs")
	

def build_phylogeny(target_name, alignment_file, fast_boot_iters=1000):
	"""
	
	Bootstrap support method information:
	http://www.iqtree.org/doc/Frequently-Asked-Questions#how-do-i-interpret-ultrafast-bootstrap-ufboot-support-values
	"""
	
	
	run("iqtree "
		f"-s {alignment_file} "
		f"-bb {fast_boot_iters} "
		f">> {console_log_name} 2>&1",
		shell=True)
	
	
	log_result(f"\n{target_name} Bootstraps", fast_boot_iters)
	with open(f"{alignment_file}.iqtree") as tree_results:
		for line in tree_results:
			if line.startswith("Random seed number"):
				log_result(f"{target_name} seed", line[20:-1])			
			elif line.startswith("Best-fit model according to BIC"):
				log_result(f"{target_name} model", line[33:-1])

	

	#Clean up additional IQtree output. I suspect essentially all users are just going to want to see a final tree,
	#so I'm tossing these files for now. In future probably should be exposed via a "--keep_intermediates" flag.
	#Also, in testing I rarely see difference b/w consensus tree and best tree. Going with best for now.
	exts_to_remove = ('*.ckp.gz', '*.iqtree', '*.contree', '*.splits.nex', '*.bionj', '*.mldist', '*.model.gz', '*.aln.log')
	for ext in exts_to_remove:
		for unused_file in glob.glob(ext):
			os.remove(unused_file)



if __name__ == "__main__":

	log_result("NTM Validation Pipeline", None, "header")
	log_result("Sample name", args.sample_name)
	log_result("Date", run_name.split('_')[-2])
	log_result("Time", run_name.split('_')[-1].replace('-',':'))
	log_result("Version", __version__)



	log_result("Read QC and filtering with fastp", None, "header")
	log_result("fastp version", 'fastp -v', 'version')
	trimmed_reads = fastp_filtering(args.nanopore_reads, args.min_fastq_length, args.min_fastq_qscore)

	log_result("De novo assembley with Flye", None, "header")
	log_result("Flye version", 'flye -v', 'version')
	flye_assembly(trimmed_reads)

	log_result("Species ANI identification with fastANI", None, "header")
	log_result("fastANI version", 'fastANI -v', 'version')
	fastani_compute(f"{run_name}_assembly.fasta")

	log_result("Target sequence extraction via BLASTn", None, "header")
	log_result("BLASTn version", 'blastn -version', 'version')
	extract_target(f"{run_name}_assembly.fasta", args.target_seq_dir, args.NTM_ref_dir, qcov_hsp_perc=args.blast_coverage)

	log_result("Target sequence alignment with Muscle", None, "header")
	log_result("Muscle version", 'muscle -version', 'version')
	for target_seq in Path(args.target_seq_dir).glob("*.fasta"):
		target_name = target_seq.stem
		ref_target_dir = f"{args.NTM_ref_dir}/{target_name}/"
		concatenate_align(f"{run_name}_assembly_{target_name}.fasta", ref_target_dir)

	log_result("ML tree phylogeny with IQ-TREE", None, "header")
	log_result("IQ-TREE version", 'iqtree -version', 'version')
	for target_seq in Path(args.target_seq_dir).glob("*.fasta"):
		target_name = target_seq.stem
		build_phylogeny(target_name, f"{run_name}_assembly_{target_name}_references.aln")


	# log_result("Tree plotting with ", None, "header")
	# log_result("IQ-TREE version", 'iqtree -version', 'version')
	for target_seq in Path(args.target_seq_dir).glob("*.fasta"):
		target_name = target_seq.stem
		newick_tree_plotter.plot_tree(
			f"{run_name}_assembly_{target_name}_references.aln.treefile",
			f"{run_name}_assembly_{target_name}.pdf",
			sample_name=f"{args.sample_name}_{target_name}",
			plot_title=target_name)





"""
Todo:

Multigene model in IQtree:
http://www.iqtree.org/doc/Advanced-Tutorial





"""











