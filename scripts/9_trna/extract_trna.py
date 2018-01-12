#!/usr/bin/python

# Script number:				9.1
# File:							1 of 4
# Prerequisite script(s):
# Prerequisite file(s):
# Description:					Extract trna from trna file downloaded from tRNADB-CE database (http://trna.ie.niigata-u.ac.jp/cgi-bin/trnadb/index.cgi)
# Output files:

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
# VARIABLES
##########################

codon_list = []
for base1 in sorted(['A', 'C', 'G', 'T']):
	for base2 in sorted(['A', 'C', 'G', 'T']):
		for base3 in sorted(['A', 'C', 'G', 'T']):
			codon_list.append(base1+base2+base3)



##########################
# FUNCTIONS
##########################


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




def extract_trna():

	trna = {}

	for seq_record in SeqIO.parse("raw_files/trna/trna_sequence.fasta", "fasta"):
		splits_description = seq_record.description.split('|')
		accession = splits_description[1]

		if accession not in trna:
			trna[accession] = []

		trna_record = [seq_record.description, seq_record.seq]
		trna[accession].append(trna_record)

	return(trna)

def write_trna_files(accession, trna):


	initiator = False
	elongator = False
	output_line = ''

	for entry in trna:
		anticodon = entry[0].split('|')[7]
		if anticodon == 'CAT':
			if 'initiator' in entry[0]:
				initiator = True
			if 'elongator' in entry[0]:
				elongator = True

		output_line += '>{}\n{}\n'.format(entry[0], entry[1])

	if initiator and elongator:
		output_file = open('outputs/trna/trna_extractions/{}_trna.fasta'.format(accession), 'w')
		output_file.write(output_line)
		output_file.close()


def write_to_files(accessions, acc_counts, trna):

	for accession in accessions:
		print('%s: %s' % (acc_counts[accession], accession))

		if accession in trna:
			write_trna_files(accession, trna[accession])


def run(accessions, acc_counts):

	trna = extract_trna()
	write_to_files(accessions, acc_counts, trna)


def main():

	t0_main = time.time()

	setupStrictDirectories('outputs/trna/')
	setupStrictDirectories('outputs/trna/trna_extractions/')

	# Return the accessions of the genomes in the directory
	accession_list = get_accessions('outputs/genome_extractions/')

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


	t1_main = time.time()
	print ('\n%s\nTotal time: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))


if __name__ == '__main__':
	main()
