#!/usr/bin/python

# Script number:				na
# File:							3 of 4
# Prerequisite script(s):
# Prerequisite file(s):			get_genomes, write_genome_files
# Description:				    Run CodonW
# Output files:


import numpy as np
import re
import time
import multiprocessing as mp
import os
import sys
import itertools
import pandas as pd
import shutil

##########################
# VARIABLES
##########################

seed_files = ["cai.coa", "cbi.coa", "coa_raw", "codon.coa", "cusort.coa", "eigen.coa", "fop.coa", "genes.coa", "hilo.coa", "summary.coa"]
output_extensions = ['.blk', '.out']

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
		setupStrictDirectory(directory)

# Get a list of the accessions
def get_folder_accessions(directory):

    accessions = []
    for folder in os.listdir(directory):
        if not folder.endswith('.DS_Store'):
            accessions.append(folder)

    return accessions

def move_files(files, directory_to):

	for file in files:
		shutil.move(file, directory_to)



def run_genomes(accessions, acc_counts):

	for accession in accessions:

		print('%s: %s' % (acc_counts[accession], accession))

		directory = 'outputs/cai//accessions/{}/'
		seed_directory = 'outputs/cai//accessions/{}/seed_outputs/'.format(accession)
		high_exp_directory = 'outputs/cai//accessions/{}/high_outputs/'.format(accession)
		exp_directory = 'outputs/cai//accessions/{}/outputs/'.format(accession)
		setupStrictDirectories([seed_directory, high_exp_directory, exp_directory])

		exp_files = ['outputs/cai//accessions/{}/{}_calc_cai_set{}'.format(accession,accession, extension) for extension in output_extensions]
		high_exp_files = ['outputs/cai//accessions/{}/{}_calc_cai_set_high{}'.format(accession,accession, extension) for extension in output_extensions]

		# Input files
		acc_file = 'outputs/cai//accessions/{}/{}_calc_cai_set.txt'.format(accession, accession)
		high_acc_file = 'outputs/cai//accessions/{}/{}_calc_cai_set_high.txt'.format(accession, accession)
		input_list = {
			'cds_file': acc_file,
			'fop': '{}/fop.coa'.format(seed_directory),
			'cai': '{}/cai.coa'.format(seed_directory),
			'cbi': '{}/cbi.coa'.format(seed_directory),
		}

		# Run seed files
		run_command1 = "packages/codonW/./codonw %s -coa_cu -coa_num %s -nomenu -silent" % (high_acc_file, "100%")
		os.system(run_command1)

		move_files(seed_files, seed_directory)
		move_files(high_exp_files, high_exp_directory)

		run_command2 = "packages/codonW/./codonw %s -all_indices -fop_file %s -cai_file %s -cbi_file %s -nomenu -silent" % (input_list['cds_file'], input_list['fop'], input_list['cai'], input_list['cbi'])
		os.system(run_command2)

		move_files(exp_files, exp_directory)



def run(accessions, acc_counts):

	# run_in_parralell(accessions, [acc_counts], run_genomes)
	run_genomes(accessions, acc_counts)



def main():


	accession_list = get_folder_accessions('outputs/cai/accessions/')

	accessions = []
	acc_counts = {}

	accession_count = 0

	for accession in sorted(accession_list):
		accession_count += 1
		if accession_count:
		# if accession_count <= 1:
			accessions.append(accession)
			acc_counts[accession] = accession_count

	run(accessions, acc_counts)


if __name__ == "__main__":
	main()
