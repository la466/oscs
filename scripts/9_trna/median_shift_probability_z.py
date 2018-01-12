#!/usr/bin/python

# Script number:				9.4
# File:							4 of 4
# Prerequisite script(s):		median_shift_probability, synonymous_site_simulation
# Prerequisite file(s):
# Description:					Group the shift probabilities with Z scores
# Output files:					median_shift_probability_z.csv

import numpy as np
import re
import time
import multiprocessing as mp
import os
import sys
import itertools
import shutil
from Bio import SeqIO
from scipy.stats import spearmanr, kendalltau
import unittest

##########################
# VARIABLES
##########################

codon_map = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
	"TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
	"TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
	"TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
	"CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
	"CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
	"CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
	"CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
	"ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
	"ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
	"AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
	"AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
	"GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
	"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
	"GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
	"GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

codon_map_t4 = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
	"TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
	"TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
	"TGT":"C", "TGC":"C", "TGA":"W", "TGG":"W",
	"CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
	"CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
	"CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
	"CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
	"ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
	"ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
	"AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
	"AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
	"GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
	"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
	"GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
	"GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

aa_map = {
	"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
	"Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
	"Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
	"Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"
}

codon_list = [codon for codon in sorted(codon_map)]
aa_list = list(set([codon_map[codon] for codon in codon_map]))
aa_list_t4 = list(set([codon_map_t4[codon] for codon in codon_map_t4]))

aa_codon_map = {}
for aa in aa_list:
	aa_codon_map[aa] = []
for codon in codon_map:
	aa_codon_map[codon_map[codon]].append(codon)


wc_pairs = {
	"A": "T",
	"C": "G",
	"G": "C",
	"T": "A"
}

wobble_pairs = {
	"A": ["C", "A", "G"],
	"C": [],
	"G": ["T"],
	"T": ["G", "T", "C"]
}

stable_wobble_pairs = {
	"A": ["C"],
	"C": [],
	"G": ["T"],
	"T": ["G", "T"]
}

unstable_wobble_pairs = {
	"A": ["A", "G"],
	"C": [],
	"G": [],
	"T": ["C"]
}


stops = ['TAA', 'TAG', 'TGA']

def get_t4_genomes():

	t4_genomes = []
	with open('outputs/gene_filtering/table_4_genomes_in_list.txt', 'rU') as t4_file:
		lines = t4_file.readlines()
		for line in lines:
			acc = line.strip('\n')
			t4_genomes.append(acc)
	return(t4_genomes)

t4_genomes = get_t4_genomes()

##########################
# FUNCTIONS
##########################

# Parrallelise
def run_in_parralell(input_list, args, function_to_run, required_workers = None):

	# blank list to hold results
	results = []

	"""
	Setup multiprocessing
	"""

	# Determine the number of cpus to use if the number of workers isnt defined
	if not required_workers:
		required_workers = (mp.cpu_count()) - 1

	# Split the list you wish to iterate over into smaller chunks that can be parallelised
	chunked_input_list = [input_list[i::required_workers] for i in range(required_workers)]

	# Multiprocessing setup
	pool = mp.Pool(required_workers)


	"""
	Iterate function over each of the chunked lists
	"""

	for i in chunked_input_list:
		current_args = args.copy()
		new_args = [i]

		for arg in current_args:
			new_args.append(arg)

		process = pool.apply_async(function_to_run, new_args)
		results.append(process)

	"""
	Close the multiprocessing threads
	"""

	pool.close()
	pool.join()


	return(results)



def get_genome_stops(accession):

	if accession in t4_genomes:
		stops = ['TAA', 'TAG']
	else:
		stops = ['TAA', 'TAG', 'TGA']

	return(stops)


# Set up the new directories for the good files and the bad files
def setupDirectories(direc):

	# If the directory already exists, delete and make new directory
	if not os.path.exists(direc):
		os.makedirs(direc)

def setupStrictDirectory(directory):
	if os.path.exists(directory):
		shutil.rmtree(directory)
	os.makedirs(directory)

def setupStrictDirectories(directory_list):

	if isinstance(directory_list, list):
		for directory in directory_list:
			setupStrictDirectory(directory)
	else:
		setupStrictDirectory(directory_list)


def get_z_scores():

	p1_zs = {
		'cs': {},
		'ss': {},
	}
	m1_zs = {
		'cs': {},
		'ss': {},
	}

	with open('outputs/simulation_synonymous_site_analysis/combined_stops.csv', 'rU') as my_file:

		lines = my_file.readlines()

		for line in lines[1:]:

			splits = line.split(',')
			acc = splits[0]
			p1_z = splits[3]
			m1_z = splits[6]

			p1_zs['ss'][acc] = p1_z
			m1_zs['ss'][acc] = m1_z

	with open('outputs/simulation_codon_shuffle_analysis/combined_stops.csv', 'rU') as my_file:

		lines = my_file.readlines()

		for line in lines[1:]:

			splits = line.split(',')
			acc = splits[0]
			p1_z = splits[3]
			m1_z = splits[6]

			p1_zs['cs'][acc] = p1_z
			m1_zs['cs'][acc] = m1_z

	return(p1_zs, m1_zs)

def get_median_costs():

	median_prob_list_p1 = {}
	median_prob_list_m1 = {}
	gcs = {}

	with open('outputs/trna/median_shift_probability.csv', 'rU') as my_file:

		lines = my_file.readlines()

		for line in lines[1:]:

			line = line.strip('\n')
			splits = line.split(',')
			acc = splits[0]
			gc = splits[1]
			gc3 = splits[2]

			gcs[acc] = [gc, gc3]

			if acc not in median_prob_list_p1:
				median_prob_list_p1[acc] = []
				median_prob_list_m1[acc] = []

			p1_p = splits[3]
			m1_p = splits[4]

			median_prob_list_p1[acc].append(float(p1_p))
			median_prob_list_m1[acc].append(float(m1_p))


	return(median_prob_list_p1, median_prob_list_m1, gcs)


def run_extracts():

	outputs = []

	p1_zs, m1_zs = get_z_scores()
	median_prob_list_p1, median_prob_list_m1, gcs = get_median_costs()

	output = [p1_zs, m1_zs, median_prob_list_p1, median_prob_list_m1, gcs]
	outputs.append(output)

	return(outputs)

def prepare_results(results, parallel, output_file_list):

	output_list = []

	if parallel:
		for result in results:
			outputs = result.get()
			for output in outputs:
				output_list.append(output)
	else:
		for output in results:
			output_list.append(output)

	output_lines = {}
	for file in output_file_list:
		output_lines[file] = []


	for output in output_list:

		p1_zs = output[0]
		m1_zs = output[1]
		median_prob_list_p1 = output[2]
		median_prob_list_m1 = output[3]
		gcs = output[4]

		for accession in sorted(median_prob_list_p1):

			# output_line2 = '{},{},{},{},{},{},{}\n'.format(accession, p1_zs['cs'][accession], p1_zs['ss'][accession], np.median(median_prob_list_p1[accession]), m1_zs['ss'][accession], m1_zs['ss'][accession], np.median(median_prob_list_m1[accession]))
			# output_lines[1].append(output_line2)

			for i in range(len(median_prob_list_p1[accession])):
			# print(median_prob_list_p1[accession])
				output_line1 = '{},{},{},{},{},{},{},{},{}\n'.format(accession, gcs[accession][0], gcs[accession][1], p1_zs['cs'][accession], p1_zs['ss'][accession], median_prob_list_p1[accession][i], m1_zs['cs'][accession], m1_zs['ss'][accession], median_prob_list_m1[accession][i])
				output_lines[0].append(output_line1)

	return(output_lines)


def write_to_files(output_file_list, output_lines):

	for file in output_file_list:

		output_file = open('outputs/trna/{}'.format(output_file_list[file][0]), 'w')
		output_file.write(output_file_list[file][1])

		for line in output_lines[file]:
			output_file.write(line)
		output_file.close()



def run():

	output_file_list = {
		0: ['median_shift_probability_z.csv', 'acc,gc,gc3,z_p1_cs,z_p1_ss,prob_p1,z_m1_cs,z_m1_ss,prob_m1\n'],
		# 1: ['median_costs_median_z.csv', 'acc,z_p1,prob_p1,z_m1,prob_m1\n'],
	}

	results = run_extracts()
	output_lines = prepare_results(results, False, output_file_list)

	# results = run_in_parralell(accessions, [acc_counts], run_genomes)
	# output_lines = prepare_results(results, True, output_file_list)

	write_to_files(output_file_list, output_lines)

def main():

	t0_main = time.time()

	run()

	t1_main = time.time()
	print ('\n%s\nTotal time: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))






if __name__ == '__main__':
	main()
