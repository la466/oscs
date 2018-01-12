#!/usr/bin/python

# Script number:				9.2
# File:							2 of 4
# Prerequisite script(s):		extract_trna
# Prerequisite file(s):
# Description:					Calculate the correlation between the costs of frameshift in CDS and OSC density for each genome
# Output files:					cds_cost_osc_density_correlations.csv

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

aa_map = {
	"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
	"Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
	"Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
	"Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"
}

codon_list = [codon for codon in sorted(codon_map)]
aa_list = list(set([codon_map[codon] for codon in codon_map]))

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


def get_anticodons(accession):

	anticodons = []
	anticodon_full_info = []
	record_id = False

	for seq_record in SeqIO.parse('outputs/trna/trna_extractions/{}_trna.fasta'.format(accession), "fasta"):

		if not record_id:
			record_id = seq_record.id

		record_description = seq_record.description
		description_splits = record_description.split('|')

		identifier = description_splits[0]
		trna_aa = description_splits[6]
		anticodon = description_splits[7]
		additional_information = description_splits[14]

		anticodons.append(anticodon)
		anticodon_full_info.append([anticodon, trna_aa, additional_information])

	unique_anticodons = list(set(anticodons))

	return(anticodons, unique_anticodons, anticodon_full_info)


def tidy_anticodon_with_description(anticodon, description):

	# Do some tidying of the methionine anticodons
	if description:
		if 'initiator' in description:
			anticodon += '|i'
		elif 'elongator' in description:
			anticodon += '|e'
		elif 'converted' in description:
			anticodon += '|{}'.format('tRNA-Met(Ile)')

	return(anticodon)


def map_aa_anticodons(anticodon_full_info):

	acc_aa_anticodon_map = {}
	for aa in aa_list:
		acc_aa_anticodon_map[aa] = []

	# For each anticodon entry inthe list
	for entry in anticodon_full_info:

		anticodon = entry[0]
		amino_acid = entry[1]
		description = entry[2]

		if description:
			anticodon = tidy_anticodon_with_description(anticodon, description)

		# Convert amino acid to letter, add to dictionary
		if amino_acid in aa_map:
			acc_aa_anticodon_map[aa_map[amino_acid]].append(anticodon)

	return(acc_aa_anticodon_map)


def check_all_amino_acids_have_anticodon(aa_anticodon_map):

	all_have_anticodons = True

	for aa in sorted(aa_anticodon_map):

		# Ignore the stop codons
		if '*' != aa:
			if all_have_anticodons and len(aa_anticodon_map[aa]) < 1:
				all_have_anticodons = False

	return(all_have_anticodons)



def check_wc_pair(nt1, nt2):

	if nt2 == wc_pairs[nt1]:
		return(True)
	else:
		return(False)


def clean_anticodon(anticodon):

	return(anticodon.split('|')[0])




def check_cognate(anticodon, codon, aa):

	# Check to see whether the codon carries the correct amino acid  corresonding
	# to the anticodon amino acid

	if aa == codon_map[codon]:

		anticodon = clean_anticodon(anticodon)

		# Reverse anticodon for nt matching
		reverse = anticodon[::-1]

		# Check to see if the first two nts for wc pairs
		pos1 = check_wc_pair(reverse[0], codon[0])
		pos2 = check_wc_pair(reverse[1], codon[1])

		# Check if the wobble pair is a wobble pair or a wc pair
		if codon[2] in wobble_pairs[reverse[2]] or codon[2] == wc_pairs[reverse[2]]:
			pos3 = True
		else:
			pos3 = False

		# If the first 2 bases form watson-crick pairs, and the third forms a wobble or wc-pair, return true
		if pos1 and pos2 and pos3:
			return(True)
		else:
			return(False)

	else:

		# Doesnt carry the correct amino acid
		return(False)


def check_cognate_stability(anticodon, codon):

	anticodon = clean_anticodon(anticodon)

	# Reverse anticodon for nt matching
	reverse = anticodon[::-1]

	# if the pairing is one of the unstable pairings
	if codon[2] in stable_wobble_pairs[reverse[2]] or codon[2] == wc_pairs[reverse[2]]:
		return(True)
	else:
		return(False)



def check_near_cognate(anticodon, codon):

	anticodon = clean_anticodon(anticodon)

	# Reverse anticodon for nt matching
	reverse = anticodon[::-1]

	pos1 = check_wc_pair(reverse[0], codon[0])
	pos2 = check_wc_pair(reverse[1], codon[1])

	if pos1 and pos2:
		return(True)

	elif pos1:
		if codon[2] in stable_wobble_pairs[reverse[2]] or codon[2] in wc_pairs[reverse[2]]:
			return(True)
		else:
			return(False)

	elif pos2:
		if codon[2] in stable_wobble_pairs[reverse[2]] or codon[2] in wc_pairs[reverse[2]]:
			return(True)
		else:
			return(False)

	else:
		return(False)




def get_codon_bindings(anticodon_set, genome_stops):

	anticodon = anticodon_set[0]
	aa = anticodon_set[1]

	binding = {
		'cognate_stable': [],
		'cognate_unstable': [],
		'near_cognate': []
	}

	splits = anticodon.split('|')

	# Special case of modified CAU
	if len(splits) > 1 and 'tRNA-Met(Ile)' in splits[1]:
		binding['cognate_stable'].append('ATA')

	else:

		for codon in codon_map:
			if codon not in genome_stops:

				# Check if carrying the correct amino acid
				cognate = check_cognate(anticodon, codon, aa) ##

				if cognate:
					# Get the stability of the wobble pair
					stable = check_cognate_stability(anticodon, codon) ##

					if stable:
						binding['cognate_stable'].append(codon)
					else:
						binding['cognate_unstable'].append(codon)

				else:
					near_cognate = check_near_cognate(anticodon, codon) ##
					if near_cognate:
						binding['near_cognate'].append(codon)

	# Sort for unittests
	for stability in binding:
		binding[stability] = sorted(binding[stability])

	return(binding)


def codon_anticodons(bindings):

	codon_anticodon_strategy = {}
	for codon in codon_list:
		codon_anticodon_strategy[codon] = {
			'cognate_stable': [],
			'cognate_unstable': [],
			'near_cognate': []
		}

	for anticodon in sorted(bindings):
		for anticodon_binding_class in bindings[anticodon]:
			codon_set = bindings[anticodon][anticodon_binding_class]

			for codon in codon_set:
				codon_anticodon_strategy[codon][anticodon_binding_class].append(anticodon)

	return(codon_anticodon_strategy)


def anticodon_strategy(acc_aa_anticodon_map, genome_stops):

	anticodon_codon_map = {}

	for aa in sorted(acc_aa_anticodon_map):
		anticodons = list(set(acc_aa_anticodon_map[aa]))

		for anticodon in sorted(anticodons):

			anticodon_splits = anticodon.split('|')

			# Only consider regular anticodons and the elongator CAT anticodon
			if len(anticodon_splits) < 2 or 'i' != anticodon_splits[1]:
				binding = get_codon_bindings([anticodon, aa], genome_stops)
				anticodon_codon_map[anticodon] = binding

	return(anticodon_codon_map)



def clean_anticodons(anticodon_full_info):

	anticodons_clean = []

	for entry in anticodon_full_info:
		anticodon = entry[0]
		if entry[2]:
			if 'initiator' in entry[2]:
				anticodon += '|i'
			elif 'elongator' in entry[2]:
				anticodon += '|e'
			elif 'converted' in entry[2]:
				anticodon += '|{}'.format('tRNA-Met(Ile)')

		if '|i' not in anticodon:
			anticodons_clean.append(anticodon)

	return(anticodons_clean)

def get_anticodon_bindings(accession, genome_stops):

	# Get a list of anticodons, unique anticodons, anticodons plus additional information
	anticodons, unique_anticodons, anticodon_full_info = get_anticodons(accession)

	# Return a map linking amino acids to the anticodons that decdode them
	aa_anticodon_map = map_aa_anticodons(anticodon_full_info)

	# Check if all amino acids have at least one anticodon
	all_amino_acids_have_anticodons = check_all_amino_acids_have_anticodon(aa_anticodon_map) ##

	# Check if methionines are present
	if not all_amino_acids_have_anticodons:
		print('Not a full set of AA\'s can be described')
		return(None)

	else:

		# Work out the cognate / near cognate rules for the anticodon set
		anticodon_codon_map = anticodon_strategy(aa_anticodon_map, genome_stops)

		# Organise into codon: anticodon set
		codon_anticodon_strategy = codon_anticodons(anticodon_codon_map)

		anticodons_clean = clean_anticodons(anticodon_full_info)

		return([codon_anticodon_strategy, anticodons_clean])



def check_bindings_output_for_cognate(bindings_output):

	codon_anticodon_strategy = bindings_output[0]

	cognate_check = True

	# Check to see if there is a cognate anticodon for every codon
	for codon in codon_anticodon_strategy:
		cognate_count = len(codon_anticodon_strategy[codon]['cognate_stable']) + len(codon_anticodon_strategy[codon]['cognate_unstable'])

		if codon not in stops and cognate_check and cognate_count == 0:
			cognate_check = False

	return(cognate_check)


def anticodon_copy_count(anticodons_clean, query_list):

	anticodon_count = 0

	for query in query_list:
		for anticodon in query:
			anticodon_count += anticodons_clean.count(anticodon)

	return(anticodon_count)



def get_codon_p(codon_anticodon_strategy, anticodons_clean, focal_codon, fs_codon):

	codon_cognate_stable = codon_anticodon_strategy[focal_codon]['cognate_stable']
	codon_cognate_unstable = codon_anticodon_strategy[focal_codon]['cognate_unstable']
	codon_near_cognates = codon_anticodon_strategy[focal_codon]['near_cognate']

	fs_cognate_stable = codon_anticodon_strategy[fs_codon]['cognate_stable']
	fs_cognate_unstable = codon_anticodon_strategy[fs_codon]['cognate_unstable']
	fs_near_cognates = codon_anticodon_strategy[fs_codon]['near_cognate']

	# Get the number of cognate anticodons for the focal codon
	cognate_count = anticodon_copy_count(anticodons_clean, [codon_cognate_stable, codon_cognate_unstable]) ##

	# Get the sets of unstable cognate or near cognate anticodons that could slip or those that cant slip
	slippery_set, unslippery_set = get_slippery_set(codon_cognate_unstable, codon_near_cognates, fs_cognate_stable, fs_cognate_unstable, fs_near_cognates) ##

	# Get the anticodon copy count for each set
	slippery_set_count = anticodon_copy_count(anticodons_clean, [slippery_set]) ##
	unslippery_set_count = anticodon_copy_count(anticodons_clean, [unslippery_set]) ##

	fs_p = calc_fs_p(cognate_count, slippery_set_count, unslippery_set_count) ##

	return(fs_p)


def get_slippery_set(focal_unstable, focal_near, shifted_stable, shifted_unstable, shifted_near):

	slippery_set = []
	unslippery_set = []

	for anticodon in focal_unstable:
		if anticodon in shifted_stable or anticodon in shifted_unstable or anticodon in shifted_near:
			slippery_set.append(anticodon)
		else:
			unslippery_set.append(anticodon)

	for anticodon in focal_near:
		if anticodon in shifted_stable or anticodon in shifted_unstable or anticodon in shifted_near:
			slippery_set.append(anticodon)
		else:
			unslippery_set.append(anticodon)


	return(slippery_set, unslippery_set)


def calc_fs_p(cognate_count, slippery_set_count, unslippery_set_count):

	b = 0.01

	numerator = b*slippery_set_count
	denominator = (b*slippery_set_count) + (b*unslippery_set_count) + cognate_count

	p = np.divide(numerator, denominator)

	return(p)

def calc_osc_proportions(cds_seqs):

	osc_proportions = {}

	genome_p1_count = 0
	genome_m1_count = 0
	genome_p1_codons = 0
	genome_m1_codons = 0

	i = 0
	for cds in cds_seqs:
		i+=1
		osc_proportions[i] = {
			'p1': [],
			'm1': [],
			'total': []
		}

		p1 = re.findall(r'.{3}', cds[4:])
		m1 = re.findall(r'.{3}', cds[2:])

		p1_count = 0
		m1_count = 0

		for stop in stops:
			p1_count += p1.count(stop)
			m1_count += m1.count(stop)

			genome_p1_count += p1.count(stop)
			genome_m1_count += m1.count(stop)

		# total = p1_count + m1_count

		p1_prop = np.divide(100, len(p1))*p1_count
		m1_prop = np.divide(100, len(m1))*m1_count
		total_prop = p1_prop + m1_prop

		osc_proportions[i]['p1'].append(p1_prop)
		osc_proportions[i]['m1'].append(m1_prop)
		# osc_proportions[i]['total'].append(total_prop)

		genome_p1_codons += len(p1)
		genome_m1_codons += len(m1)


	genome_total = genome_p1_count + genome_m1_count
	genome_p1_density = (100/genome_p1_codons)*genome_p1_count
	genome_m1_density = (100/genome_m1_codons)*genome_m1_count
	# genome_total_density = genome_p1_density + genome_m1_density

	genome_osc_densities = {}
	genome_osc_densities['p1'] = genome_p1_density
	genome_osc_densities['m1'] = genome_m1_density
	# genome_osc_densities['total'] = genome_total_density

	return(genome_osc_densities)


def calc_seq_ps(cds, codon_anticodon_strategy, anticodons_clean):

	ps = []
	for i in range(3, len(cds)-3, 3):
		focal_codon = cds[i:i+3]
		fsp1_codon = cds[i+1:i+4]
		p = get_codon_p(codon_anticodon_strategy, anticodons_clean, focal_codon, fsp1_codon)
		ps.append(p)

	return(ps)



def calc_genome_total_p(ps, genome_osc_densities):

	cors = {}

	cds_ps = []

	for cds in ps:
		cds_ps.append(sum(ps[cds]))

	genome_p = sum(cds_ps)


	return(genome_p)


def calc_post_p1(i, cds):

	query_sequence = cds[i:]
	p1_codons = re.findall(r'.{3}', query_sequence[1:])
	found_stop = False

	count = 0
	for codon in p1_codons:
		if not found_stop:
			if codon not in stops:
				count +=1
			else:
				found_stop = True


	return(count)

def calc_post_m1(i, cds):

	query_sequence = cds[i:]
	p1_codons = re.findall(r'.{3}', query_sequence[2:])
	found_stop = False

	count = 0
	for codon in p1_codons:
		if not found_stop:
			if codon not in stops:
				count +=1
			else:
				found_stop = True


	return(count)


def calc_nprenpost(cds):

	p1posts = []
	m1posts = []
	p1preposts = []
	m1preposts = []


	for i in range(3, len(cds)-3, 3):
		pre = int(i/3)

		p1post = calc_post_p1(i, cds)
		p1prepost = pre+p1post

		m1post = calc_post_m1(i-3, cds)
		m1prepost = pre+m1post

		p1posts.append(p1post)
		p1preposts.append(p1prepost)
		m1posts.append(m1post)
		m1preposts.append(m1prepost)



	return(p1posts, p1preposts, m1posts, m1preposts)


def calc_codon_costs(ps, preposts):

	codon_costs = {"p1": [], "m1": []}

	p1_preposts = preposts[2]
	m1_preposts = preposts[3]

	for i in range(len(ps)):
		codon_costs["p1"].append(ps[i]*p1_preposts[i])
		codon_costs["m1"].append(ps[i]*m1_preposts[i])

	return(codon_costs)


def calc_osc_densities(cds):

	p1_count = 0
	m1_count = 0

	p1_codons = re.findall(r'.{3}', cds[4:])
	m1_codons = re.findall(r'.{3}', cds[2:])

	for stop in stops:
		p1_count += p1_codons.count(stop)
		m1_count += m1_codons.count(stop)

	p1_density = (100/len(p1_codons))*p1_count
	m1_density = (100/len(m1_codons))*m1_count

	return(p1_density, m1_density)


def sequence_analysis(accession, cds_seqs, bindings_output):

	codon_anticodon_strategy = bindings_output[0]
	anticodons_clean = bindings_output[1]

	genome_osc_densities = calc_osc_proportions(cds_seqs[:])

	p1_densities = []
	m1_densities = []
	p1_costs = []
	m1_costs = []

	for cds in cds_seqs[:]:
		ps = calc_seq_ps(cds, codon_anticodon_strategy, anticodons_clean)
		preposts = calc_nprenpost(cds)
		codon_costs = calc_codon_costs(ps, preposts)
		p1_density, m1_density = calc_osc_densities(cds)

		p1_costs.append(sum(codon_costs["p1"]))
		p1_densities.append(p1_density)
		m1_costs.append(sum(codon_costs["m1"]))
		m1_densities.append(m1_density)

	corp1 = spearmanr(p1_costs, p1_densities)
	corm1 = spearmanr(m1_costs, m1_densities)

	return(corp1, corm1, genome_osc_densities)



def run_genomes(accessions, acc_counts):

	outputs = []

	for accession in accessions:

		print('%s: %s' % (acc_counts[accession], accession))

		cds_seqs = get_seqs(accession) ##
		gc, gc3 = calc_gc(cds_seqs) ##

		genome_stops = get_genome_stops(accession)

		bindings_output = get_anticodon_bindings(accession, genome_stops) ##

		# If the bindings output exists
		if bindings_output:

			cognate_check = check_bindings_output_for_cognate(bindings_output) ##

			if cognate_check:
				corp1, corm1, genome_osc_densities = sequence_analysis(accession, cds_seqs, bindings_output)
				# cors_max_p_per_cds, genome_p, genome_osc_densities = sequence_analysis(accession, cds_seqs, bindings_output)
				output = [accession, gc, gc3, corp1, corm1, genome_osc_densities]
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
		accession = output[0]
		gc = output[1]
		gc3 = output[2]
		corp1 = output[3]
		corm1 = output[4]
		genome_osc_densities = output[5]

		output_line1 = '{},{},{},{},{},{},{}\n'.format(accession, gc, gc3, corp1.correlation, corp1.pvalue, corm1.correlation, corm1.pvalue)
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

	accession_list = get_accessions('outputs/trna/trna_extractions/')

	accessions = []
	acc_counts = {}

	accession_count = 0

	for accession in sorted(accession_list):
		accession_count += 1
		if accession_count:
		# if accession == 'AE017263':
		# if accession_count <= 1:
			accessions.append(accession)
			acc_counts[accession] = accession_count


	# results = run_genomes(accessions, acc_counts)
	# output_lines = prepare_results(results, False)

	results = run_in_parralell(accessions, [acc_counts], run_genomes)


	output_file_list = {
		0: ['cds_cost_osc_density_correlations.csv', 'acc,gc,gc3,p1_cor,p1_p,m1_cor,m1_p\n'],
	}

	output_lines = prepare_results(results, True, output_file_list)
	write_to_files(output_file_list, output_lines)

def main():

	t0_main = time.time()

	run()

	t1_main = time.time()
	print ('\n%s\nTotal time: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))






if __name__ == '__main__':
	main()
