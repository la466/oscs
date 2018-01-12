#!/usr/bin/python

# Script number:				2.1
# File:							1 of 1
# Prerequisite script(s):
# Prerequisite file(s):
# Description:					Sorts genomes

import os
import sys
import imp
import time
import re
import shutil
from datetime import timedelta
from Bio import SeqIO


#############
# Variables #
#############

script_name = "2.1_sort_genomes.py\n"
script_description = "Sort genomes, leaving only 1 per genus > 500,000 BP\n"

embl_directory = "raw_files/bacteria_raw/"


avoid = ['unverified', 'uncultured', 'alpha', 'beta', 'gamma']

# Avoided accessions with no cds information
avoided = ['AP014723', 'AP014507', 'AP014638']

output_directory = "outputs/"
output_genome_sorting = output_directory + "genome_sorting/"
output_good_genomes = output_genome_sorting + "good_genomes/"
output_t4_genomes = output_genome_sorting + "t4_genomes/"

output_genome_list_path = output_genome_sorting + 'genome_list.csv'
output_genome_list_path_t4 = output_genome_sorting + 'genome_list_t4.csv'


#############
# Functions #
#############

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


def get_files(directory, filetype):

	files = []

	# For each file in the directory
	for file in os.listdir(directory):

		# If the file ends with the filetype
		if file.endswith(filetype):

			# Append the file to the files list
			files.append(file)

	return files, len(files)

# Clean the name of the genus
def clean_genus_name(genus_name):

	characters = ["[", "]", "'"]

	for character in characters:
		if character in genus_name:
			genus_name = genus_name.replace(character, '')

	return genus_name


def read_genome(genome_path, acc, genomes, t4_genomes):

	with open(genome_path, 'rU') as genome_file:

		description = ''
		bp = 0
		OS = 0
		ID = 0
		CDS = 0

		transl_table = False
		got_transl_table = False

		for line in genome_file:

			if line.startswith('OS') and OS == 0:

				description = line.replace('OS   ', '')		# Remove the 'OS' tag from the line
				description = description.replace('\n', '')		# Remove the line break
				OS += 1

			if line.startswith('ID') and ID == 0:
				genome_bp = re.findall('(?<=; )\d+(?= BP)', line)
				bp = int(genome_bp[0])
				ID += 1

			if line.startswith('FT   CDS'):
				CDS += 1

			if not got_transl_table:
				if line.startswith('FT'):
					if 'transl_table' in line:
						transl_table = re.findall(r'transl_table=(\d+)', line)[0]
						got_transl_table = True


		description_split = description.split(' ')
		genus = description_split[0]

		candidatus = False
		if description_split[0] == 'Candidatus':
			genus = description_split[1]
			candidatus = True
		else:
			genus = description_split[0]


		genus = clean_genus_name(genus)		# Tidy genus string

		if genus not in genomes and bp > 500000:			# If genus has not already been selected
			if genus not in avoid and genus[0].isupper() and CDS > 0:	# Check if genome should be avoided and is true
				if candidatus == True:
					genus = 'Candidatus %s' % genus
				genomes[genus] = [bp, description, acc]


		if transl_table == '4' and bp > 500000:			# If genus has not already been selected
			if genus not in avoid and genus[0].isupper() and CDS > 0:	# Check if genome should be avoided and is true
				if candidatus == True:
					genus = 'Candidatus %s' % genus
				t4_genomes[acc] = [acc, genus, description, bp]



def run():

	start_time = time.time()	# Set the start time

	script_misc() 	# Description

	# Create directories
	create_directory(output_directory)
	create_directory(output_genome_sorting)
	create_strict_directory(output_good_genomes)
	create_directory(output_t4_genomes)

	raw_embl_accessions, number_raw_embl_accessions = get_files(embl_directory, ".embl")	# Get a list of the raw EMBL files

	genome_count = 0
	genomes = {}
	t4_genomes = {}

	# Examine the genomes and get accessions of unique genuses
	for raw_acc in raw_embl_accessions:

		genome_count += 1
		acc = raw_acc.strip(".embl")

		# if genome_count > 1257 and genome_count < 1260:
		print (acc)
		print ('Genome %s of %s' % (genome_count, number_raw_embl_accessions))
		print ('-' * 20)

		if acc not in avoided:
			genome_path = embl_directory + raw_acc
			read_genome(genome_path, acc, genomes, t4_genomes)

	print ('Copying genomes')

	# Output genome info to file and copy to good genomes folder
	output_genome_list = open(output_genome_list_path, 'w')
	output_genome_list.write('genus,description,base_pairs,accession\n')

	for genome in genomes:
		output_string = "%s,%s,%s,%s\n" % (genome, genomes[genome][1], genomes[genome][0], genomes[genome][2])
		output_genome_list.write(output_string)

		print ('Copying %s to %s' % (genomes[genome][2], output_good_genomes))
		acc_path = embl_directory + genomes[genome][2] + '.embl'
		shutil.copy(acc_path, output_good_genomes)

	output_genome_list_t4 = open(output_genome_list_path_t4, 'w')
	output_genome_list_t4.write('genus,description,base_pairs,accession\n')

	for genome in t4_genomes:
		output_string = "%s,%s,%s,%s\n" % (t4_genomes[genome][1], t4_genomes[genome][2], t4_genomes[genome][3],genome)
		output_genome_list_t4.write(output_string)

		print ('Copying %s to %s' % (genome, output_t4_genomes))
		acc_path = embl_directory + genome + '.embl'
		shutil.copy(acc_path, output_t4_genomes)



	# Print the elapsed time
	elapsed = time.time() - start_time
	print (str(timedelta(seconds=elapsed)))


###############
# Run script #
###############

if __name__ == "__main__":
	run()
