#!/usr/bin/python

# Script number:				3.1
# File:							1 of 1
# Prerequisite script(s):
# Prerequisite file(s):
# Description:					Filter and write CDSs in fasta format
# Output files:					filter_errors.csv, filter_errors_t4.csv, table_4_genomes_in_list.txt

import os
import sys
import imp
import time
import re
import shutil
from datetime import timedelta
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio.Seq import Seq
import multiprocessing as mp


#############
# Variables #
#############

script_name =  os.path.basename(__file__)
script_description = "Convert embl file to fasta format. Filter genes.\n"

embl_directory = "raw_files/bacteria_raw/"
output_directory = "outputs/"

genome_sorting = output_directory + "genome_sorting/"
good_genomes = genome_sorting + "good_genomes/"
t4_genomes = genome_sorting + "t4_genomes/"

output_fasta_directory = output_directory + 'genome_extractions/'
output_fasta_directory_t4 = output_directory + 'genome_extractions_t4/'
output_filter_directory = output_directory + 'gene_filtering/'


genome_list = genome_sorting + 'genome_list.csv'
genome_list_t4 = genome_sorting + 'genome_list_t4.csv'

output_error_file_path = output_filter_directory + 'filter_errors.csv'
output_error_file_path_t4 = output_filter_directory + 'filter_errors_t4.csv'

output_table_4_path = output_filter_directory + 'table_4_genomes_in_list.txt'

#############
# Functions #
#############

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


# Create new directory
def create_directory(directory_path):

	if os.path.exists(directory_path):
		print ('Directory already exists: %s\n' % directory_path)
	else:
		print ('Making new directory: %s\n' % directory_path)
		os.mkdir(directory_path)


# Strictly create new directory (remove it previsouly)
def create_strict_directory(directory_path):

	if os.path.exists(directory_path):
		print ('Directory already exists: %s' % directory_path)
		print ('Removing directory')
		shutil.rmtree(directory_path)
		print ('Making new directory: %s\n' % directory_path)
		os.mkdir(directory_path)
	else:
		print ('Making new directory: %s\n' % directory_path)
		os.mkdir(directory_path)


# Get the script descriptions
def script_misc():
	print(script_name)
	print(script_description)


def get_accession_list(genome_list_path):

	accession_list = {}
	header = True

	with open(genome_list_path, 'rU') as genome_file:
		for entry in genome_file:
			if header:
				header = False
				continue

			entry_splits = entry.split(',')
			accession = entry_splits[3].strip('\n')
			accession_list[accession] = [entry_splits[0], entry_splits[1], entry_splits[2]]

	return accession_list


def get_sequence_from_raw(genome_sequence, cds_location, cds_start, cds_end, cds_strand, operator):

	cds_location_nice = ''
	cds_locations_nice = []

	if operator == 'join':

		cds_sequence = ''
		cds_locations_nice = []

		cds_locations = re.findall('(?<=\{).+(?=})', str(cds_location))[0]
		cds_locations = cds_locations.split(',')
		for location in cds_locations:
			location = location.strip(' ')

			cds_limits = location.split(':')
			join_start = int(re.findall('\d+', cds_limits[0])[0])
			join_end = int(re.findall('\d+', cds_limits[1])[0])
			join_strand = re.findall('(?<=\().+(?=\))', location)[0]

			join_location_nice = '%s:%s(%s);' % (join_start+1, join_end,join_strand)
			cds_locations_nice.append(join_location_nice)

			if join_strand == '-':
				cds_sequence += genome_sequence[join_start:join_end].reverse_complement()
			else:
				cds_sequence += genome_sequence[join_start:join_end]

		for join_location_nice in reversed(cds_locations_nice):
			cds_location_nice += join_location_nice

	else:
		if cds_strand == -1:
			cds_sequence = genome_sequence[cds_start:cds_end].reverse_complement()
		else:
			cds_sequence = genome_sequence[cds_start:cds_end]

		cds_location_nice += '%s:%s' % (cds_start+1,cds_end)


	return cds_sequence, cds_location_nice


def check_length_three(cds):

	if len(cds) % 3 != 0:
		filter_check = False
	else:
		filter_check = True

	return filter_check

def check_bases(cds):

	filter_check = True
	bases = ['A', 'C', 'T', 'G']
	for nt in cds:
		if nt not in bases:
			filter_check = False

	return filter_check

def check_standard_stop(cds, translation_table):

	filter_check = True
	stop_codons = {}
	stop_codons[11] = ['TAA', 'TGA', 'TAG']
	stop_codons[4] = ['TAA', 'TAG']

	if cds[-3:] not in stop_codons[int(translation_table)]:
		filter_check = False

	return filter_check


def check_inframe_stop(cds, translation_table):

	filter_check = True
	stop_codons = {}
	stop_codons[11] = ['TAA', 'TGA', 'TAG']
	stop_codons[4] = ['TAA', 'TAG']

	for i in range(0, len(cds)-3, 3):
		if cds[i:i+3] in stop_codons[int(translation_table)]:
			filter_check = False

	return filter_check


def parse_genome(accession, record_path):

	records = 0

	output_record = ''
	error_string = ''

	translation_table_4 = False


	for seq_record in SeqIO.parse(record_path, "embl"):

		genome_sequence = seq_record.seq


		for seq_feature in seq_record.features:															# For each of the genome features
			if seq_feature.type=="CDS":																	# Check if feature is a CDS

				try:
					len(seq_feature.qualifiers['translation'])==1									# Check translation exists
					records += 1

					# Extract CDS information
					protein_id = seq_feature.qualifiers['protein_id'][0]									# Protein ID
					transl_table = seq_feature.qualifiers['transl_table'][0]								# Translation table
					operator = seq_feature.location_operator												# Location operator (join)
					cds_strand = seq_feature.strand															# Strand
					cds_start = seq_feature.location.nofuzzy_start
					cds_end = seq_feature.location.nofuzzy_end
					cds_location = seq_feature.location

					if transl_table == '4':
						translation_table_4 = True

					try:
						pseudo = seq_feature.qualifiers['pseudo']
					except:
						pass


					try:
						locus_tag = seq_feature.qualifiers['locus_tag'][0]
					except:
						locus_tag = 'no_locus_tag_cds%s' % records

					try:																					# See if the cds name exists
						cds_name = seq_feature.qualifiers['gene'][0]										# Gene name
					except:
						cds_name = 'no_cds_qualifier_%s' % seq_feature.qualifiers['locus_tag'][0]

					cds_sequence, cds_location_nice = get_sequence_from_raw(genome_sequence, cds_location, cds_start, cds_end, cds_strand, operator)		# CDS sequence


					# Filter cds
					cds_multiple_three_check = check_length_three(cds_sequence)
					cds_base_check = check_bases(cds_sequence)
					cds_check_stop = check_standard_stop(cds_sequence, transl_table)
					cds_inframe_stop = check_inframe_stop(cds_sequence, transl_table)

					if not cds_multiple_three_check:
						error_string += '%s,%s,1\n' % (accession, locus_tag)
					if not cds_base_check:
						error_string += '%s,%s,2\n' % (accession, locus_tag)
					if not cds_check_stop:
						error_string += '%s,%s,3\n' % (accession, locus_tag)
					if not cds_inframe_stop:
						error_string += '%s,%s,4\n' % (accession, locus_tag)

					cds_entry = ''
					if cds_multiple_three_check and cds_base_check and cds_check_stop and cds_inframe_stop:
						cds_entry = '>%s|%s|%s|%s|%s|%s|%s\n' % (accession, locus_tag, protein_id, cds_name, transl_table, cds_location_nice, cds_strand)
						cds_entry += '%s\n' % cds_sequence

						output_record += cds_entry

				except:
					pass

	return output_record, error_string, translation_table_4

def extract_cds(accession, table):

	if table == 11:
		record_output_path = output_fasta_directory + accession + '.txt'
		record_path = good_genomes + accession + '.embl'
	elif table == 4:
		record_output_path = output_fasta_directory_t4 + accession + '.txt'
		record_path = t4_genomes + accession + '.embl'



	try:
		output_record, error_string, table4 = parse_genome(accession, record_path)
	except:

		# Workaround for BioPython thrown up by awkward annotation
		# https://github.com/biopython/biopython/issues/341

		line_count = 0
		f = open(record_path, "rU")
		lines = f.readlines()
		f.close()

		output_lines = []
		for line in lines:
		    if line.startswith('CO'):
		       	output_lines.append("#" + line)
		    else:
		    	output_lines.append(line)

		f = open(record_path, "w")
		f.write("".join(output_lines))
		f.close()

		# Now parse the genome
		output_record, error_string, table4 = parse_genome(accession, record_path)

	return output_record, record_output_path, error_string, accession, table4



def write_to_file(output_record, record_output_path):

	print ('Writing to: %s' % record_output_path)

	output_file = open(record_output_path, 'w')
	output_file.write(output_record)
	output_file.close()




def run_genomes(accessions, acc_counts, table):

	outputs = []

	for accession in accessions:

		print ('{}: {}'.format(accession, acc_counts[accession]))
		output_record, record_output_path, error_string, accession, table4 = extract_cds(accession, table)
		output = [output_record, record_output_path, error_string, accession, table4]

		outputs.append(output)

	return (outputs)




def run():

	start_time = time.time()	# Set the start time

	script_misc() 	# Description

	print('Processing good genomes')

	create_strict_directory(output_fasta_directory)
	create_strict_directory(output_filter_directory)

	accession_list = get_accession_list(genome_list)

	accessions = []
	acc_counts = {}

	accession_count = 0

	for accession in sorted(accession_list):
		accession_count += 1
		if accession_count:
		# if accession_count <= 10:
			accessions.append(accession)
			acc_counts[accession] = accession_count

	results = run_in_parralell(accessions, [acc_counts, 11], run_genomes)


	error_file = open(output_error_file_path, 'w')
	error_file.write('acc,locus_tag,error_code,,1=length_fail,2=non_actg,3=non_standard_stop,4=in_frame_stop\n')

	table_4_file = open(output_table_4_path, 'w')

	for result in results:
		outputs = result.get()
		for output in outputs:
			output_record = output[0]
			record_output_path = output[1]
			error_string = output[2]
			accession = output[3]
			table4 = output[4]

			write_to_file(output_record, record_output_path)
			error_file.write(error_string)

			if table4:
				line = '%s\n' % accession
				table_4_file.write(line)

	error_file.close()
	table_4_file.close()


	print('Processing good genomes')

	create_strict_directory(output_fasta_directory_t4)
	accession_list = get_accession_list(genome_list_t4)

	accessions = []
	acc_counts = {}

	accession_count = 0

	for accession in sorted(accession_list):
		accession_count += 1
		if accession_count:
			accessions.append(accession)
			acc_counts[accession] = accession_count

	results = run_in_parralell(accessions, [acc_counts, 4], run_genomes)

	error_file = open(output_error_file_path_t4, 'w')
	error_file.write('acc,locus_tag,error_code,,1=length_fail,2=non_actg,3=non_standard_stop,4=in_frame_stop\n')

	for result in results:
		outputs = result.get()
		for output in outputs:
			output_record = output[0]
			record_output_path = output[1]
			error_string = output[2]
			accession = output[3]
			table = output[4]

			write_to_file(output_record, record_output_path)
			error_file.write(error_string)


	error_file.close()


	# Print the elapsed time
	elapsed = time.time() - start_time
	print (str(timedelta(seconds=elapsed)))


###############
# Run script #
###############

if __name__ == "__main__":
	run()
