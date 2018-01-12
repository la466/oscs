#!/usr/bin/python

# Script number:				1.1
# File:							1 of 1
# Prerequisite script(s):
# Prerequisite file(s): bacteria_accessions.txt
# Description: Downloads genomes from EMBL using accessions provided

import os
import sys
import imp
import urllib2
import HTMLParser
import time
from datetime import timedelta

#############
# Variables #
#############

script_name = "1.1_get_genomes.py\n"
script_description = "Downloading bacteria genomes from EMBL website\n"

url = "http://www.ebi.ac.uk/Tools/dbfetch/emblfetch?db=embl&id=AE005175&format=embl&style=raw&Retrieve=Retrieve"

output_directory = "raw_files/"
bacteria_output_directory = output_directory + "bacteria_raw/"

accessions_path = "misc/bacteria_accessions.txt"


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


# Get a list of the accession numbers
def read_accessions(file_path):

	accession_list = []

	accession_file = open(file_path, 'U')
	for line in accession_file:
		line = line.replace("\n", "")
		accession_list.append(line)

	return accession_list, len(accession_list)


# Download the genomes from EMBL
def download_genomes(bacteria_accession_list, number_accessions, start_time):

	acc_count = 0

	for acc in bacteria_accession_list:

		acc_count += 1
		print("Retrieving %s: %s of %s") % (acc, acc_count, number_accessions)

		embl_path = "http://www.ebi.ac.uk/ena/data/view/" + acc + "&display=txt&expanded=true"

		embl_response = urllib2.urlopen(embl_path)
		embl_data = embl_response.read()

		# Write to file
		output_file_path = bacteria_output_directory + acc + ".embl"
		output_file = open(output_file_path, 'w')
		output_file.write(embl_data)
		output_file.close

		# Print the elapsed time
		elapsed = time.time() - start_time
		print (str(timedelta(seconds=elapsed)))


###############
# Run script #
###############


def run():

	start_time = time.time()

	script_misc() 	# Description
	bacteria_accession_list, number_accessions = read_accessions(accessions_path)	# Get a list of the accession numbers

	# Create the output directories
	create_strict_directory(output_directory)
	create_strict_directory(bacteria_output_directory)

	download_genomes(bacteria_accession_list, number_accessions, start_time)	# Download bacteria EMBL files



if __name__ == "__main__":
	run()
