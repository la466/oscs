#!/usr/bin/python

# Script number:				6.5
# File:							5 of 6
# Prerequisite script(s):       53mm
# Prerequisite file(s):
# Description:					Analyse Markov model simulations
# Output files:                 stops.csv, all_codons.csv, combined_stops.csv, combined_stops_percentage_excess.csv


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

# Map of codons to their amino acids
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
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

# Sorted list of codons
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

# List of stop codons
stop_codons = ['TAA', 'TAG', 'TGA']

# Stop codons defined by translation table
stops = {}
stops[4] = ['TAA', 'TAG']
stops[11] = ['TAA', 'TAG', 'TGA']

# Reading frames
frames = [1,2,'both']


def get_t4_genomes():

	t4_genomes = []
	with open('outputs/gene_filtering/table_4_genomes_in_list.txt', 'rU') as t4_file:
		lines = t4_file.readlines()
		for line in lines:
			acc = line.strip('\n')
			t4_genomes.append(acc)
	return(t4_genomes)

t4 = get_t4_genomes()

output_directory = 'outputs/mmodels/53mm/'

##########################
# FUNCTIONS
##########################


def create_directory(directory_path):

	if os.path.exists(directory_path):
		print ('Directory already exists: %s\n' % directory_path)
	else:
		print ('Making new directory: %s\n' % directory_path)
		os.mkdir(directory_path)


def get_file_paths(directory):

    file_paths = {}
    for file in os.listdir(directory):
        if file.endswith('.csv'):
            file_paths[re.findall(r'[^_]*', file)[0]] = file

    return file_paths


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


# Functions
def get_stops(accession):

    if accession in t4:
        genome_stops = stops[4]
    else:
        genome_stops = stops[11]

    return(genome_stops)


def get_lines(file_path):

    acc_codons = []
    real_counts = []
    randomised_counts = []

    with open(output_directory + file_path, 'rU') as randomised_file:

        lines = randomised_file.readlines()

        codon_count = 0

        line_count = 0
        for line in lines:
            line_count += 1

            if line_count == 1:
                line = line.strip('\n')
                splits = line.split(',')
                for i in range(4, 67):
                    acc_codons.append(splits[i][:-2])


            if line_count == 2:
                real_counts.append(line.strip('\n').split(','))
                gc = float(line.split(',')[1])
                gc3 = float(line.split(',')[2])
                codon_count = float(line.split(',')[3])
            elif line_count > 2:
                randomised_counts.append(line.strip('\n').split(','))

    return(acc_codons, real_counts, randomised_counts, gc, gc3, codon_count)



def counts_to_dict(accession, acc_codons, real_counts, randomised_counts, codon_count):

    real_counts_dict = {}
    randomised_counts_dict = {}

    for frame in frames:
        real_counts_dict[frame] = {}
        randomised_counts_dict[frame] = {}
        for codon in codon_list:
            real_counts_dict[frame][codon] = 0
            randomised_counts_dict[frame][codon] = []

    for i in range(len(codon_list)):

        real_counts_dict[1][codon_list[i]] = float(real_counts[0][i+4])
        real_counts_dict[2][codon_list[i]] = float(real_counts[0][i+68])
        real_counts_dict['both'][codon_list[i]] = float(real_counts[0][i+4])+ float(real_counts[0][i+68])


    for repeat in randomised_counts:
        for i in range(len(codon_list)):
            randomised_counts_dict[1][codon_list[i]].append(float(repeat[i+4]))
            randomised_counts_dict[2][codon_list[i]].append(float(repeat[i+68]))
            randomised_counts_dict['both'][codon_list[i]].append(float(repeat[i+4])+float(repeat[i+68]))

    return(real_counts_dict, randomised_counts_dict)


def analyse_counts(acc_codons, real_counts_dict, randomised_counts_dict, genome_stops):

    z_scores = {}
    pvals = {}

    for frame in frames:
        z_scores[frame] = {}
        pvals[frame] = {}

    for frame in real_counts_dict:
        for codon in codon_list:

            z = (real_counts_dict[frame][codon] - np.mean(randomised_counts_dict[frame][codon])) / np.std(randomised_counts_dict[frame][codon])
            z_scores[frame][codon] = z
            pvals[frame][codon] = scipy.stats.ttest_1samp(randomised_counts_dict[frame][codon], real_counts_dict[frame][codon])



    # Combine the stop codons
    real_totals = {}
    randomised_totals = {}
    for frame in real_counts_dict:
        total_real = 0
        for stop in genome_stops:
            total_real += real_counts_dict[frame][stop]
        real_totals[frame] = total_real

        randomised_totals[frame] = []
        for i in range(len(randomised_counts_dict[frame]['TAA'])):
            total_randomised = 0
            for stop in genome_stops:
                total_randomised += randomised_counts_dict[frame][stop][i]
            randomised_totals[frame].append(total_randomised)

    z_totals = {}
    pes = {}

    for frame in real_totals:
        z = (real_totals[frame] - np.mean(randomised_totals[frame])) / np.std(randomised_totals[frame])

        ttest = scipy.stats.ttest_1samp(randomised_totals[frame], real_totals[frame])
        z_totals[frame] = [z,ttest]

        pes[frame] = (real_totals[frame] - np.mean(randomised_totals[frame]))/real_totals[frame]

    return(z_scores, pvals, z_totals, pes)


# Test function to check columns are being returned correctly
def check_values(real_counts):

    check = True

    if float("{0:.9f}".format(real_counts[1]['AAA'])) != 6.348513179:
        check = False
        print('AAA 1 fail')
    if float("{0:.9f}".format(real_counts[2]['AAA'])) != 4.821302629:
        check = False
        print('AAA 2 fail')
    if float("{0:.9f}".format(real_counts[1]['TTT'])) != 4.101674622:
        check = False
        print('TTT 1 fail')
    if float("{0:.9f}".format(real_counts[2]['TTT'])) != 6.967885319:
        check = False
        print('TTT 2 fail')

    if not check:
        print('Column checks failed')

    return(check)



def run_genome(accessions, file_paths, acc_counts):

    outputs = []

    for accession in accessions:

        t0 = time.time()

        print('%s: %s' % (acc_counts[accession], accession))

        # Get genome stop codons
        genome_stops = get_stops(accession)


        # Return the lines from the file
        acc_codons, real_counts, randomised_counts, gc, gc3, codon_count = get_lines(file_paths[accession])

        # Group the counts in a dictionary
        real_counts_dict, randomised_counts_dict = counts_to_dict(accession, acc_codons, real_counts, randomised_counts, codon_count)

        # Analyse the counts
        z_scores, pvals, z_totals, pes = analyse_counts(acc_codons, real_counts_dict, randomised_counts_dict, genome_stops)
        # codon_z_scores, codon_le_scores, combined_stops_z_scores, combined_stops_le_scores, codon_t_tests, combined_stops_t_tests = analyse_counts(acc_codons, real_counts_dict, randomised_counts_dict, genome_stops)

        if accession == 'AE000511':
            check = check_values(real_counts_dict)
            if not check:
                break

        output = [accession, gc, gc3, genome_stops, z_scores, pvals, z_totals, pes]
        outputs.append(output)

        t1 = time.time()
        # print('Genome time: %s' % (t1-t0))

    return (outputs)



def write_stop_codons(results):

    output_file = open('outputs/mmodels_analysis/53mm/stops.csv', 'w')
    header = 'acc,gc,gc3'
    for frame in frames:
        for stop in sorted(stop_codons):
            header += ',%s_%s_z,%s_%s_pval' % (stop, frame, stop, frame)
    header += '\n'
    output_file.write(header)


    for result in results:

        outputs = result.get()
        for output in outputs:

            acc = output[0]
            gc = output[1]
            gc3 = output[2]
            genome_stops = output[3]
            z_scores = output[4]
            pvals = output[5]

            output_line = '%s,%s,%s' % (acc, gc, gc3)
            for frame in z_scores:
                for stop in sorted(stop_codons):
                    output_line += ',%s,%s' % (z_scores[frame][stop], pvals[frame][stop].pvalue)
            output_line += '\n'
            output_file.write(output_line)

    output_file.close()

def write_combined_stops(results):

    output_file = open('outputs/mmodels_analysis/53mm/combined_stops.csv', 'w')
    header = 'acc,gc,gc3'
    for frame in frames:
        header += ',osc_%s_z,osc_%s_ttest,osc_%s_pval' % (frame, frame,frame)
    header += '\n'
    output_file.write(header)


    for result in results:

        outputs = result.get()
        for output in outputs:

            acc = output[0]
            gc = output[1]
            gc3 = output[2]
            combined_zs = output[6]


            output_line = '%s,%s,%s' % (acc, gc, gc3)
            for frame in frames:
                output_line += ',%s,%s,%s' % (combined_zs[frame][0], combined_zs[frame][1].statistic, combined_zs[frame][1].pvalue)

            output_line += '\n'
            output_file.write(output_line)

    output_file.close()


def write_all_codons(results):

    output_file = open('outputs/mmodels_analysis/53mm/all_codons.csv', 'w')
    header = 'acc,gc,gc3'
    for frame in frames:
        for codon in sorted(codon_list):
            header += ',%s_%s_z,%s_%s_pval' % (codon, frame, codon, frame)
    header += '\n'
    output_file.write(header)


    for result in results:

        outputs = result.get()
        for output in outputs:

            acc = output[0]
            gc = output[1]
            gc3 = output[2]
            genome_stops = output[3]
            z_scores = output[4]
            pvals = output[5]

            output_line = '%s,%s,%s' % (acc, gc, gc3)
            for frame in z_scores:
                for codon in sorted(codon_list):
                    output_line += ',%s,%s' % (z_scores[frame][codon], pvals[frame][codon].pvalue)
            output_line += '\n'
            output_file.write(output_line)

    output_file.close()

def write_percentage_excess(results):

    output_file = open('outputs/mmodels_analysis/53mm/combined_stops_percentage_excess.csv', 'w')
    header = 'acc,gc,gc3'
    for frame in frames:
        header += ',pe_%s' % (frame)
    header += '\n'
    output_file.write(header)


    for result in results:

        outputs = result.get()
        for output in outputs:

            acc = output[0]
            gc = output[1]
            gc3 = output[2]
            pes = output[7]

            output_line = '%s,%s,%s' % (acc, gc, gc3)
            for frame in pes:
                output_line += ',%s' % (pes[frame])
            output_line += '\n'
            output_file.write(output_line)

    output_file.close()


def main():

    t0_main = time.time()

    create_directory('outputs/mmodels_analysis/')
    create_directory('outputs/mmodels_analysis/53mm/')


    file_paths = get_file_paths(output_directory)

    accessions = []
    acc_counts = {}

    accession_count = 0
    for accession in sorted(file_paths):
        accession_count += 1
        if accession_count:
        # if accession_count <= 1:
            accessions.append(accession)
            acc_counts[accession] = accession_count


    # results = run_genome(accessions, file_paths, acc_counts)
    results = run_in_parralell(accessions, [file_paths, acc_counts], run_genome)

    write_stop_codons(results)
    write_combined_stops(results)
    write_all_codons(results)
    write_percentage_excess(results)




    t1_main = time.time()
    print ('\n%s\nTotal time: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))


if __name__ == '__main__':
    main()
