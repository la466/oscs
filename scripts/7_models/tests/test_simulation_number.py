#!/usr/bin/python

# Script number:			    7 Extra2
# File:
# Prerequisite script(s):       ecoli_codon_shuffle_500_repeats
# Prerequisite file(s):
# Description:					Get differernt repeat numbers of codon shuffle for E coli
# Output files:                 means.csv

import numpy as np
import re
import time
import multiprocessing as mp
import os
import sys
from scipy.stats import t
import scipy.stats

##########################
# Variables
##########################

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

codon_list = [codon for codon in sorted(codon_map)]

# Get a sorted list of codons for each amino acid
aa_map = {}
for codon in codon_map:
    if codon_map[codon] not in aa_map:
        aa_map[codon_map[codon]] = [codon]
    else:
        aa_map[codon_map[codon]].append(codon)
for aa in aa_map:
    aa_map[aa].sort()


stop_codons = ['TAA', 'TAG', 'TGA']
stops = {}
stops[4] = ['TAA', 'TAG']
stops[11] = ['TAA', 'TAG', 'TGA']

frames = [1,2,'both']

tgn_codons = ['TAA', 'TAC', 'TAG', 'TAT']

def get_t4_genomes():

	t4_genomes = []
	with open('outputs/gene_filtering/table_4_genomes_in_list.txt', 'rU') as t4_file:
		lines = t4_file.readlines()
		for line in lines:
			acc = line.strip('\n')
			t4_genomes.append(acc)
	return(t4_genomes)

t4 = get_t4_genomes()

##########################
# FUNCTIONS
##########################


def create_directory(directory_path):

	if os.path.exists(directory_path):
		print ('Directory already exists: %s\n' % directory_path)
	else:
		print ('Making new directory: %s\n' % directory_path)
		os.mkdir(directory_path)

def get_accessions(directory):

    accessions = []
    for file in os.listdir(directory):
        if not file.endswith('.DS_Store'):
            accessions.append(re.findall(r'(.+)(?=_shuffle_codons\.csv)', file)[0])

    return accessions



def run_in_parralell(input_list, args, function_to_run, workers = None, onebyone = False):

    if not workers:
        workers = (mp.cpu_count()) - 2

    if not onebyone:
        chunk_list = [input_list[i::workers] for i in range(workers)]
    else:
        chunk_list = input_list

    pool = mp.Pool(workers)
    results = []

    for i in chunk_list:
        current_args = args.copy()
        new_args = [i]



        for arg in current_args:
            # print(arg)
            new_args.append(arg)

        process = pool.apply_async(function_to_run, new_args)
        results.append(process)

    pool.close()
    pool.join()

    return(results)


def get_stops(accession):

    if accession in t4:
        genome_stops = stops[4]
    else:
        genome_stops = stops[11]

    return(genome_stops)


def get_lines(accession):

    acc_codons = []
    real_counts = []
    randomised_counts = []

    with open('outputs/simulation_tests/{}_codon_shuffle.csv'.format(accession), 'rU') as randomise_file:

        lines = randomise_file.readlines()

        line_count = 0
        for line in lines:
            line_count += 1

            if line_count == 1:
                line = line.strip('\n')
                splits = line.split(',')
                for i in range(4, 68):
                    acc_codons.append(splits[i][:-2])


            if line_count == 2:
                real_counts.append(line.strip('\n').split(','))
                gc = float(line.split(',')[1])
                gc3 = float(line.split(',')[2])
                codon_count = float(line.split(',')[3])
            elif line_count > 2:
                randomised_counts.append(line.strip('\n').split(','))


    return(acc_codons, real_counts, randomised_counts, gc, gc3, float(codon_count))



def counts_to_dict(acc_codons, real_counts, randomised_counts, codon_count):

    real_counts_dict = {}
    randomised_counts_dict = {}

    for frame in frames:
        real_counts_dict[frame] = {}
        randomised_counts_dict[frame] = {}
        for codon in acc_codons:
            real_counts_dict[frame][codon] = 0
            randomised_counts_dict[frame][codon] = []


    for i in range(len(acc_codons)):
        # print(i, real_counts[0][i+4])
        real_counts_dict[1][acc_codons[i]] = (100/codon_count)*float(real_counts[0][i+4])
        real_counts_dict[2][acc_codons[i]] = (100/codon_count)*float(real_counts[0][i+68])
        real_counts_dict['both'][acc_codons[i]] = (100/codon_count)*(float(real_counts[0][i+4])+ float(real_counts[0][i+68]))


    for repeat in randomised_counts:
        for i in range(len(acc_codons)):
            randomised_counts_dict[1][acc_codons[i]].append((100/codon_count)*float(repeat[i+4]))
            randomised_counts_dict[2][acc_codons[i]].append((100/codon_count)*float(repeat[i+68]))
            randomised_counts_dict['both'][acc_codons[i]].append((100/codon_count)*(float(repeat[i+4])+float(repeat[i+68])))


    return(real_counts_dict, randomised_counts_dict)



def analyse_randomisations(randomised_counts_dict):


    required = [5,10,25,50,75,100,150,200,300,400,500]
    stop_codons = ['TAA', 'TAG', 'TGA']
    frames = [1,2,'both']

    means = {}
    sds = {}
    for frame in frames:
        means[frame] = {}
        sds[frame] = {}
        for require in required:
            means[frame][require] = []
            sds[frame][require] = []


    # For each sample size
    for require in required:

        # Set up blank list to hold OSC count means
        output_means = []

        # 100 times
        for i in range(0, 100):

            np.random.seed()

            random_indicies = np.random.randint(0, 499, size=require)

            for frame in frames:

                osc_counts = []

                for index in random_indicies:
                    oscs = 0
                    for stop in stop_codons:
                        oscs += randomised_counts_dict[frame][stop][index]
                    osc_counts.append(oscs)

                mean = np.mean(osc_counts)
                sd = np.std(osc_counts)

                means[frame][require].append(mean)
                sds[frame][require].append(sd)



    output_file = open('outputs/simulation_tests/means.csv', 'w')

    header = 'sample'
    for frame in means:
        for require in sorted(means[frame]):
            header += ',frame_%s_%s' % (frame, require)
    header += '\n'
    output_file.write(header)


    for i in range(0,100):
        output_line = '%s' % (i+1)
        for frame in means:
            for require in sorted(means[frame]):
                output_line += ',%s' % (means[frame][require][i])

        output_line += '\n'
        output_file.write(output_line)

    output_file.close()


def test_levels():

    accession = 'AE005174'

    # Return the lines from the file
    acc_codons, real_counts, randomised_counts, gc, gc3, codon_count = get_lines(accession)

    # Group the counts in a dictionary
    real_counts_dict, randomised_counts_dict = counts_to_dict(acc_codons, real_counts, randomised_counts, codon_count)

    analyse_randomisations(randomised_counts_dict)

    t1 = time.time()





def main():

    t0_main = time.time()

    create_directory('outputs/simulation_tests/')
    test_levels()

    t1_main = time.time()
    print ('\n%s\nTotal time: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))


if __name__ == '__main__':
    main()
