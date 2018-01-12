#!/usr/bin/python

# Script number:				na
# File:							2 of 4
# Prerequisite script(s):
# Prerequisite file(s):			get_genomes
# Description:					Write genome files for use with CodonW
# Output files:				


import numpy as np
import re
import time
import multiprocessing as mp
import os
import sys
import itertools
import shutil

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


# Get a list of the accessions
def get_accessions(directory):

    accessions = []
    for file in os.listdir(directory):
        if not file.endswith('.DS_Store'):
            accessions.append(re.findall(r'(.+)(?=\.)', file)[0])

    return accessions


# Get a list of the cds sequences for a genome
def get_seqs(accession):

	cds_seqs = {}

	file_path = 'outputs/genome_extractions/%s.txt' % accession
	with open(file_path, 'rU') as my_file:

		line = my_file.read()
		cds_list = re.findall('>.+\n.+\n', line)

		for cds in cds_list:
			entry = cds.split('\n')

			# cds_seqs.append(cds[1])
			# cds_info.append(cds[0])

			cds_identifier = extract_identifier(entry[0])
			cds_seqs[cds_identifier] =  entry[1]

	return (cds_seqs)



def extract_identifier(cds_info):

	splits = cds_info.split('|')
	cds_identifier = splits[3]

	return(cds_identifier)


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


def get_half_cds_index(cds):

	# Remove the first 30 nts
	cds = cds[30:]
	codons = int(len(cds)/3)
	first_half_index = round(codons/2)*3

	return(cds, first_half_index)


def write_to_files(accession, cds_seqs, high_exp_identifiers):

	output_directory = 'outputs/cai/accessions/{}'.format(accession)

	setupDirectories(output_directory)

	output_file_calc_cai_set = open(output_directory + '/{}_calc_cai_set.txt'.format(accession), 'w')
	output_file_calc_cai_set_high = open(output_directory + '/{}_calc_cai_set_high.txt'.format(accession), 'w')
	output_file_analysis_set = open(output_directory + '/{}_analysis_set.txt'.format(accession), 'w')
	output_file_analysis_set_high = open(output_directory + '/{}_analysis_set_high.txt'.format(accession), 'w')


	for cds in cds_seqs:
		output_line = '>{}|{}\t\t{} residues\n'.format(cds, accession, len(cds_seqs[cds]))

		cds_seq, first_half_index = get_half_cds_index(cds_seqs[cds])
		cds_first_half = cds_seq[:first_half_index]
		cds_last_half = cds_seq[first_half_index:]

		if cds in high_exp_identifiers[accession]:
			# HE, CAI set
			output_line_calc_cai_set_high = output_line + '{}\n'.format(cds_first_half)
			output_file_calc_cai_set_high.write(output_line_calc_cai_set_high)

			# HE, analysis set
			output_line_analysis_set_high = output_line + '{}\n'.format(cds_last_half)
			output_file_analysis_set_high.write(output_line_analysis_set_high)
		else:
			output_line_calc_cai_set = output_line + '{}\n'.format(cds_first_half)
			output_file_calc_cai_set.write(output_line_calc_cai_set)

			output_line_analysis_set = output_line + '{}\n'.format(cds_last_half)
			output_file_analysis_set.write(output_line_analysis_set)

	output_file_calc_cai_set.close()
	output_file_calc_cai_set_high.close()
	output_file_analysis_set.close()
	output_file_analysis_set_high.close()



def run_genomes(accessions, acc_counts, high_exp_identifiers):

	for accession in accessions:

		if accession in high_exp_identifiers:

			print('%s: %s' % (acc_counts[accession], accession))

			cds_seqs = get_seqs(accession)
			write_to_files(accession, cds_seqs, high_exp_identifiers)


def run(accessions, acc_counts):

	setupDirectories('outputs/cai/')
	setupStrictDirectories('outputs/cai/accessions/')


	high_exp_identifiers = get_high_exp_identifiers()

	results = run_in_parralell(accessions, [acc_counts, high_exp_identifiers], run_genomes)
	# results = run_genomes(accessions, acc_counts, high_exp_identifiers)



def main():

    t0_main = time.time()

    # Return the accessions of the genomes in the directory
    accession_list = get_accessions('outputs/genome_extractions/')

    accessions = []
    acc_counts = {}

    accession_count = 0

    for accession in sorted(accession_list):
        accession_count += 1
        if accession_count:
        # if accession_count <= 4:
            accessions.append(accession)
            acc_counts[accession] = accession_count


    run(accessions, acc_counts)

    t1_main = time.time()
    print ('\n%s\nTotal time: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))


if __name__ == '__main__':
    main()
