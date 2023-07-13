#!/usr/bin/env python3

"""
https://www.ncbi.nlm.nih.gov/gene/886354
https://evolution.genetics.washington.edu/phylip/doc/consense.html
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC84584/

"""

#Imported modules
import os
import sys
from Bio import SeqIO
from pathlib import Path
from subprocess import check_output

#Authorship information
__author__ = "Angus Angermeyer"
__copyright__ = "Copyright 2023"
__license__ = "GPL"
__version__ = "0.1.0"
__email__ = "angus.angermeyer@gmail.com"


ref_file = sys.argv[1] # designed for multi fasta query file
genome_directory = sys.argv[2]


# ref_record = SeqIO.read(ref_file, "fasta")
# ref_length = len(ref_record)
ref_cutoff = 90

outdir = Path(ref_file).stem+"_extracted"

os.mkdir(outdir)

for genome_file in Path(genome_directory).glob("*.fasta"):
	genome_stem = genome_file.stem
	print(genome_stem)


	blast_call = ("blastn "
	f"-subject {genome_file} "
	f"-query {ref_file} "
	f"-outfmt 6 "
	f"-qcov_hsp_perc {ref_cutoff} "
	)

	blast_call = blast_call.split()
	
	blast_output = check_output(blast_call, text=True)

	hit_list = []

#	query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
	blast_output = blast_output.split("\n")
	for entry in blast_output:
		entry = entry.split("\t")
		
		hit_list.append([entry[-1],entry])
	
	hit_list.sort()	
	
	top_hit = hit_list[-1][1]
# 	print(hit_list[-1][1])


	if len(top_hit) > 1:
		print(top_hit)
		
		for seq_record in SeqIO.parse(genome_file, "fasta"):
			if seq_record.id == top_hit[1]:

				
				sstrt = int(top_hit[8])
				sstop = int(top_hit[9])
				
				if sstrt < sstop:
					sequence = str(seq_record.seq[sstrt-1:sstop]) #Weird slicing accounts for indexing differential b/w blast and biopython

					with open(f"{outdir}/{genome_stem}_{top_hit[0]}.fasta", 'a') as fasta_output:
						fasta_output.write(f">{top_hit[1]}_{top_hit[0]}\n")
						fasta_output.write(sequence+"\n")


				elif sstrt > sstop:
					sequence = str(seq_record.seq[sstop-1:sstrt].reverse_complement())

					with open(f"{outdir}/{genome_stem}_{top_hit[0]}.fasta", 'a') as fasta_output:
						fasta_output.write(f">{top_hit[1]}_{top_hit[0]}\n")
						fasta_output.write(sequence+"\n")
						
						
				print(sequence)











