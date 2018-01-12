#!/usr/bin/python

# Script number:				9.x
# File:							na
# Prerequisite script(s):		extract_trna
# Prerequisite file(s):
# Description:					Sanity check comparing with Warnceke (2010). Get the number of anticodons, anticodon gene copy number for each genome
# Output files:					anticodon_repetoire.csv


import numpy as np
import re
import time
import multiprocessing as mp
import os
import sys
import itertools
import shutil
from Bio import SeqIO


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
			accessions.append(re.findall(r'[^_, ^\.]*', file)[0])

	return accessions


def get_seqs(accession):

    file_path = 'outputs/genome_extractions/%s.txt' % accession
    with open(file_path, 'rU') as my_file:

        cds_seqs = []

        line = my_file.read()
        cds_list = re.findall('>.+\n.+\n', line)

        for cds in cds_list:
            cds = cds.split('\n')

            cds_seqs.append(cds[1])

    return cds_seqs


def calc_gc(cds_list):

    gc_count = gc3_count = nt_count = nt3_count = 0

    cds_count = 0
    codon_count = 0

    for cds in cds_list:
        cds_count += 1

        nt_count += len(cds[:-3])

        test_cds = cds
        gc_count += re.subn(r'[GC]', '0', test_cds[:-3])[1]

        pos3s = re.findall(r'.{2}(.)', cds[3:-3])
        gc3_count += (pos3s.count('G') + pos3s.count('C'))

        codon_count += len(pos3s)

    gc = gc_count / float(nt_count)
    gc3 = gc3_count / float(codon_count)

    return(gc,gc3)



def run_trna_file(accessions, acc_counts):

	outputs = []

	for accession in accessions:
		print('%s: %s' % (acc_counts[accession], accession))

		cds_seqs = get_seqs(accession)
		gc, gc3 = calc_gc(cds_seqs)

		anticodons = []
		record_id = False

		for seq_record in SeqIO.parse('outputs/trna/trna_extractions/{}_trna.fasta'.format(accession), "fasta"):

			if not record_id:
				record_id = seq_record.id


			record_description = seq_record.description
			description_splits = record_description.split('|')

			identifier = description_splits[0]
			trna_aa = description_splits[6]
			anticodon = description_splits[7]

			anticodons.append(anticodon)

		unique_anticodons = list(set(anticodons))

		output = [accession, gc, gc3, record_id, anticodons, unique_anticodons]
		outputs.append(output)

	return(outputs)


def list_to_string(list_of_items, delimiter):
	return('{}'.format(delimiter).join(list_of_items))

# def write_to_file_linear(outputs):
#
# 	output_file1 = open('outputs/trna/accession_anticodons.txt', 'w')
#
# 	output_file2 = open('outputs/trna/anticodon_repetoire.csv', 'w')
# 	output_file2.write('acc,gc,gc3,total_anticodons,unique_anticodons\n')
#
# 	for output in outputs:
# 		accession = output[0]
# 		gc = output[1]
# 		gc3 = output[2]
# 		record_id = output[3]
# 		anticodons = output[4]
# 		unique_anticodons = output[5]
#
# 		record_id = record_id.split('|')[1:]
#
# 		output_file1.write('>{}\n{}\n'.format(list_to_string(record_id, '|'), list_to_string(anticodons, ',')))
# 		output_file2.write('{},{},{},{},{}\n'.format(accession, gc, gc3, len(anticodons), len(unique_anticodons)))
#
#
# 	output_file1.close()
# 	output_file2.close()

def write_to_file(results):


	output_file = open('outputs/trna/anticodon_repetoire.csv', 'w')
	output_file.write('acc,gc,gc3,total_anticodons,unique_anticodons\n')

	for result in results:
		outputs = result.get()
		for output in outputs:
			accession = output[0]
			gc = output[1]
			gc3 = output[2]
			record_id = output[3]
			anticodons = output[4]
			unique_anticodons = output[5]

			record_id = record_id.split('|')[1:]

			output_file.write('{},{},{},{},{}\n'.format(accession, gc, gc3, len(anticodons), len(unique_anticodons)))


	output_file.close()



def run(accessions, acc_counts):

	results = run_in_parralell(accessions, [acc_counts], run_trna_file)
	write_to_file(results)


def main():

	t0_main = time.time()

	# Return the accessions of the genomes in the directory
	accession_list = get_accessions('outputs/trna/trna_extractions/')

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
