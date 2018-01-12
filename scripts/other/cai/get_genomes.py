#!/usr/bin/python

# Script number:				na
# File:							1 of 4
# Prerequisite script(s):
# Prerequisite file(s):
# Description:					Get genomes that contain the correct genes
# Output files:					high_exp_identifiers.txt


import numpy as np
import re
import time
import multiprocessing as mp
import os
import sys
import itertools

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


# Get a list of the accessions
def get_accessions(directory):

    accessions = []
    for file in os.listdir(directory):
        if not file.endswith('.DS_Store'):
            accessions.append(re.findall(r'(.+)(?=\.)', file)[0])

    return accessions


# Get a list of the cds sequences for a genome
def get_seqs(accession):

	file_path = 'outputs/genome_extractions/%s.txt' % accession
	with open(file_path, 'rU') as my_file:

		cds_seqs = []
		cds_info = []

		line = my_file.read()
		cds_list = re.findall('>.+\n.+\n', line)

		for cds in cds_list:
			cds = cds.split('\n')

			cds_seqs.append(cds[1])
			cds_info.append(cds[0])

	return (cds_seqs, cds_info)



def extract_identifiers(cds_info):

	cds_identifiers = []

	for cds in cds_info:
		splits = cds.split('|')
		cds_identifiers.append(splits[3])

	return(cds_identifiers)


def count_he_identifiers(cds_identifiers):

	high_identifiers = []
	for gene in highExpGenes:
		if gene in cds_identifiers and gene not in high_identifiers:
			high_identifiers.append(gene)

	return(len(high_identifiers), high_identifiers)



def run_genomes(accessions, acc_counts):

	outputs = []

	for accession in accessions:

		print('%s: %s' % (acc_counts[accession], accession))

		cds_seqs, cds_info = get_seqs(accession)
		cds_identifiers = extract_identifiers(cds_info)
		count, high_identifiers = count_he_identifiers(cds_identifiers)

		if count >= 20:
			outputs.append([accession, high_identifiers])

	return(outputs)


def return_results(results):

	setupDirectories('outputs/cai/')
	output_file = open('outputs/cai/high_exp_identifiers.txt', 'w')

	for result in results:
		outputs = result.get()
		for output in outputs:
			accession = output[0]
			high_identifiers = output[1]

			output_line = '>{}\n'.format(accession)
			output_line += '{}'.format(high_identifiers[0])
			for identifier in high_identifiers[1:20]:
				output_line += ',{}'.format(identifier)
			output_line += '\n'
			output_file.write(output_line)

	output_file.close()



def run(accessions, acc_counts):

	results = run_in_parralell(accessions, [acc_counts], run_genomes)
	# results = run_genomes(accessions, acc_counts)
	return_results(results)


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
        # if accession_count <= 10:
            accessions.append(accession)
            acc_counts[accession] = accession_count


    run(accessions, acc_counts)

    t1_main = time.time()
    print ('\n%s\nTotal time: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))


if __name__ == '__main__':
    main()
