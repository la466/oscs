#!/usr/bin/python

# Script number:				na
# File:							4 of 4
# Prerequisite script(s):
# Prerequisite file(s):			get_genomes, write_genome_files, run_codonw
# Description:				    Analyse CodonW outputs
# Output files:					cai_osc_densities.csv



import numpy as np
import re
import time
import multiprocessing as mp
import os
import sys
import itertools
import pandas as pd

##########################
# VARIABLES
##########################


genomes_dir = 'outputs/genome_extractions'

highExpGenes = [
	"rplA", "rplB", "rplC", "rplD", "rplE", "rplF", "rplI", "rplJ", "rplK", "rplL", "rplM", "rplN", "rplO", "rplP", "rplQ", "rplR", "rplS", "rplT", "rplU",
	"rpl1", "rpl2", "rpl3", "rpl4", "rpl5", "rpl6", "rpl9", "rpl10", "rpl11", "rpl12", "rpl13", "rpl14", "rpl15", "rpl16", "rpl17", "rpl18", "rpl19", "rpl20", "rpl21",
	"rpsB", "rpsC", "rpsD", "rpsE", "rpsF", "rpsG", "rpsH", "rpsI", "rpsJ", "rpsK", "rpsL", "rpsM", "rpsN", "rpsO", "rpsP", "rpsQ", "rpsR", "rpsS", "rpsT", "rpsU",
	"rps2", "rps3", "rps4", "rps5", "rps6", "rps7", "rps8", "rps9", "rps10", "rps11", "rps12", "rps13", "rps14", "rps15", "rps16", "rps17", "rps18", "rps19", "rps20", "rps21"
]

codon_map = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
	"TCT":"S", "TCC":"s", "TCA":"S", "TCG":"S",
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
	"GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}




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


# Set up the new directories for the good files and the bad files
def setupDirectories(direc):

	# If the directory already exists, delete and make new directory
	if not os.path.exists(direc):
		os.makedirs(direc)


# Get a list of the accessions
def get_accessions(directory):

	accessions = []
	for file in os.listdir(directory):
		if not file.endswith('.DS_Store'):
			accessions.append(re.findall(r'(.+)(?=\.)', file)[0])

	return accessions

# Get a list of the accessions
def get_folder_accessions(directory):

	accessions = []
	for folder in os.listdir(directory):
		if not folder.endswith('.DS_Store'):
			accessions.append(folder)

	return accessions

# Get a list of the cds sequences for a genome
def get_seqs(accession):

	file_path = 'outputs/genome_extractions/%s.txt' % accession
	with open(file_path, 'rU') as my_file:

		cds_seqs = {}

		line = my_file.read()
		cds_list = re.findall('>.+\n.+\n', line)

		for cds in cds_list:
			cds = cds.split('\n')
			splits = cds[0].split('|')
			cds_seqs[splits[3]] = cds[1]

	return (cds_seqs)

# Get a list of the cds sequences for a genome
def get_analysis_seqs(accession):

	file_path = 'outputs/cai//accessions/{}/{}_analysis_set.txt'.format(accession, accession)

	with open(file_path, 'rU') as my_file:
		cds_seqs = {}
		line = my_file.read()
		cds_list = re.findall('>.+\n.+\n', line)

		for cds in cds_list:
			cds = cds.split('\n')
			splits = cds[0].split('|')
			cds_seqs[splits[0][1:]] = cds[1]

	return (cds_seqs)

def calc_gc(cds_list):

	gc_count = gc3_count = nt_count = nt3_count = 0

	cds_count = 0
	codon_count = 0

	for cds in cds_list:
		cds_count += 1

		nt_count += len(cds_list[cds][:-3])

		test_cds = cds_list[cds]
		gc_count += re.subn(r'[GC]', '0', test_cds[:-3])[1]

		pos3s = re.findall(r'.{2}(.)', cds_list[cds][3:-3])
		gc3_count += (pos3s.count('G') + pos3s.count('C'))

		codon_count += len(pos3s)

	gc = gc_count / float(nt_count)
	gc3 = gc3_count / float(codon_count)

	return(gc,gc3)

# def extract_identifiers(cds_info):
#
#     cds_identifiers = []
#
#     for cds in cds_info:
#         splits = cds.split('|')
#         cds_identifiers.append(splits[3])
#
#     return(cds_identifiers)


def get_high_exp_identifiers():

	high_exp_identifiers = {}

	with open('outputs/cai/high_exp_identifiers.txt', 'rU') as he_file:

		lines = he_file.read()
		acc_list = re.findall('>.+\n.+\n', lines)

		for genome in acc_list:
			info = genome.split('\n')
			acc = info[0][1:]
			genes = info[1]

			high_exp_identifiers[acc] = genes

	return(high_exp_identifiers)


def extract_cais(accession, high_exp_identifiers, analysis_seqs):

	cais = {}
	ready_file = False

	with open('outputs/cai//accessions/{}/outputs/{}_calc_cai_set.out'.format(accession, accession), 'rU') as cai_file:
		lines = cai_file.readlines()

		# Arbritary number to check if CAI values have been calculated
		if len(lines) > 10:
			ready_file = True

			line_count = 0
			for line in lines:
				line_count += 1
				if line_count > 1:
					line_splits = line.split('\t')
					identifier = line_splits[0].split('|')[0]

					if 'no_cds_qualifier' in identifier:
						cai = line_splits[5]
					else:
						cai = line_splits[7]

					if identifier not in high_exp_identifiers and identifier in analysis_seqs:
						cais[identifier] = cai

	return(cais, ready_file)


def get_median_cai(cais):

	cai_list = []

	for cds in cais:
		cai_list.append(float(cds[1]))

	median = np.median(cai_list)
	return(median)


def get_osc_densities(cais, analysis_seqs):

	stop_codons = ['TAA', 'TAG', 'TGA']

	codon_count = 0
	osc_count_p1 = 0
	osc_count_m1 = 0

	for cds in cais:
		codon_count += len(analysis_seqs[cds[0]][:-3])/3

		p1_codons = re.findall(r'.{3}', analysis_seqs[cds[0]][4:])
		m1_codons = re.findall(r'.{3}', analysis_seqs[cds[0]][2:])

		for stop in stop_codons:
			osc_count_p1 += p1_codons.count(stop)
			osc_count_m1 += m1_codons.count(stop)

	osc_density_p1 = (100/codon_count) * osc_count_p1
	osc_density_m1 = (100/codon_count) * osc_count_m1

	osc_densities = {'p1': osc_density_p1, 'm1': osc_density_m1}

	return(osc_densities)




def analyse(cais, analysis_seqs):


	sorted_cais = [(k, cais[k]) for k in sorted(cais, key=cais.get, reverse=True)]
	threshold = int(round(len(sorted_cais)*0.1))

	high_cais = sorted_cais[:threshold]
	low_cais = sorted_cais[-threshold:]

	median_high_cai = get_median_cai(high_cais)
	median_low_cai = get_median_cai(low_cais)

	osc_densities_high = get_osc_densities(high_cais, analysis_seqs)
	osc_densities_low = get_osc_densities(low_cais, analysis_seqs)

	return(median_high_cai, median_low_cai, osc_densities_high, osc_densities_low)


def string_to_list(string, delimiter):
	return(string.split(delimiter))





def run_genomes(accessions, acc_counts, high_exp_identifiers):

	outputs = []

	for accession in accessions:

		print('%s: %s' % (acc_counts[accession], accession))

		accession_high_exp_identifiers = string_to_list(high_exp_identifiers[accession], ',')

		cds_seqs = get_seqs(accession)
		gc, gc3 = calc_gc(cds_seqs)

		analysis_seqs = get_analysis_seqs(accession)

		cais, read_file = extract_cais(accession, accession_high_exp_identifiers, analysis_seqs)
		if read_file:
			median_high_cai, median_low_cai, osc_densities_high, osc_densities_low = analyse(cais, analysis_seqs)
			output = [accession, gc, gc3, median_high_cai, median_low_cai, osc_densities_high, osc_densities_low]
			outputs.append(output)

	return(outputs)



def write_results(results, parallel):

	output_list = []

	if parallel:
		for result in results:
			outputs = result.get()
			for output in outputs:
				output_list.append(output)
	else:
		for output in results:
			output_list.append(output)


	output_file = open('outputs/cai/cai_osc_densities.csv', 'w')
	output_file.write('acc,gc,gc3,high_median_cai,low_median_cai,high_osc_density_p1,low_osc_density_p1,high_osc_density_m1,low_osc_density_m1\n')

	for output in output_list:
		accession = output[0]
		gc = output[1]
		gc3 = output[2]
		median_high_cai = output[3]
		median_low_cai = output[4]
		osc_densities_high = output[5]
		osc_densities_low = output[6]

		output_file.write('{},{},{},{},{},{},{},{},{}\n'.format(accession, gc, gc3, median_high_cai, median_low_cai, osc_densities_high['p1'], osc_densities_low['p1'], osc_densities_high['m1'], osc_densities_low['m1']))

	output_file.close()

def run(accessions, acc_counts, high_exp_identifiers):


	results = run_in_parralell(accessions, [acc_counts, high_exp_identifiers], run_genomes)
	write_results(results, True)

	# results = run_genomes(accessions, acc_counts, high_exp_identifiers)
	# write_results(results, False)



def main():

	t0_main = time.time()

	# Return the accessions of the genomes in the directory
	accession_list = get_folder_accessions('outputs/cai//accessions/')
	high_exp_identifiers = get_high_exp_identifiers()

	accessions = []
	acc_counts = {}

	accession_count = 0

	for accession in sorted(accession_list):
		accession_count += 1
		if accession_count:
		# if accession_count <= 1:
			accessions.append(accession)
			acc_counts[accession] = accession_count


	run(accessions, acc_counts, high_exp_identifiers)

	t1_main = time.time()
	print ('\n%s\nTotal time: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))


if __name__ == '__main__':
	main()
