#!/usr/bin/python

# Script number:				4.2
# File:							2 of 2
# Prerequisite script(s):
# Prerequisite file(s):
# Description:					Get OSC stats for table 4 genomes
# Output files:					off_frame_densities_t4.csv, combined_off_frame_densities_t4.csv

import numpy as np
import re
import time
import multiprocessing as mp
import os
import sys
import scipy.stats
import random
import collections


#############
# Variables #
#############


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

bases = ['A', 'C', 'G', 'T']

all_stops = ['TAA', 'TAG', 'TGA']

stop_codons = {}
stop_codons[4] = ['TAA', 'TAG']
stop_codons[11] = ['TAA', 'TAG', 'TGA']

t4 = ['AE015450', 'AE017263', 'AF222894', 'CP002082']

tgn_codons = ['TGA', 'TGC', 'TGG', 'TGT', 'TAA', 'TAC', 'TAG', 'TAT']

frames = [1, 2, 'both']

#############
# Functions #
#############


def run_in_parralell(input_list, args, function_to_run, workers = None, onebyone = False):

    if not workers:
        workers = (mp.cpu_count()) - 1

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


# Create new directory
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
            accessions.append(re.findall(r'(.+)(?=\.)', file)[0])

    return accessions


def get_seqs(accession):

    file_path = 'outputs/genome_extractions_t4/%s.txt' % accession
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
    for cds in cds_list:

        nts = list(cds[:-3])
        gc_count += nts.count('G')
        gc_count += nts.count('C')
        nt_count += len(nts)

        nt_index = 0
        for nt in nts:
            nt_index += 1
            if nt_index > 3 and nt_index % 3 == 0:
                if nt == 'G' or nt == 'C':
                    gc3_count += 1
                nt3_count += 1

    gc = gc_count / float(nt_count)
    gc3 = gc3_count / float(nt3_count)

    return(gc,gc3)


def get_genome_stops(accession):

    genome_stops = stop_codons[4]
    trans_table = 4


    return(genome_stops, trans_table)


def calculate_densities(cds_seqs, genome_stops):

    codon_count = 0

    # Setup blank dictionaries
    codon_counts, codon_densities = {}, {}
    for frame in frames:
        codon_counts[frame] = {}
        codon_densities[frame] = {}
        for codon in codon_list:
            codon_counts[frame][codon] = 0
            codon_densities[frame][codon] = 0

    combined_counts = collections.defaultdict(int)
    combined_densities = collections.defaultdict(int)


    for cds in cds_seqs:

        codon_count += (len(cds[:-3])/3)
        fs1 = re.findall(r'.{3}', cds[1:])
        fs2 = re.findall(r'.{3}', cds[2:])

        for codon in sorted(codon_list):
            codon_counts[1][codon] += fs1.count(codon)
            codon_counts[2][codon] += fs2.count(codon)
            codon_counts['both'][codon] += (fs1.count(codon) + fs2.count(codon))

        for stop in genome_stops:
            combined_counts[1] += fs1.count(stop)
            combined_counts[2] += fs2.count(stop)
            combined_counts['both'] += (fs1.count(stop) + fs2.count(stop))


    for frame in codon_counts:
        for codon in codon_counts[frame]:
            codon_densities[frame][codon] = (100/codon_count) * codon_counts[frame][codon]

    for frame in combined_counts:
        combined_densities[frame] = (100/codon_count) * combined_counts[frame]


    return(len(cds_seqs), codon_count, codon_densities, combined_densities)



def get_osc_densities(accession, cds_seqs):

    # Get the stop codons for the genome
    genome_stops, trans_table = get_genome_stops(accession)

    cds_count, codon_count, codon_densities, combined_densities = calculate_densities(cds_seqs, genome_stops)

    output = [cds_count, codon_count, trans_table, codon_densities, combined_densities]

    return(output)




def run_genomes(accessions, acc_counts):

    outputs = []

    for accession in accessions:

        print('%s: %s' % (acc_counts[accession], accession))

        cds_seqs = get_seqs(accession)
        gc, gc3 = calc_gc(cds_seqs)

        output = get_osc_densities(accession, cds_seqs)

        output = [accession, gc, gc3] + output

        outputs.append(output)

    return(outputs)


def write_codon_densities(results):

    output_file = open('outputs/osc_densities/off_frame_densities_t4.csv', 'w')

    header = 'acc,gc,gc3,trans_table,cds_count,codon_count'
    for frame in frames:
        for codon in sorted(codon_list):
            header += ',%s_%s' % (codon, frame)
    header += '\n'
    output_file.write(header)

    for result in results:
        outputs = result.get()
        for output in outputs:
            acc = output[0]
            gc = output[1]
            gc3 = output[2]
            cds_count = output[3]
            codon_count = output[4]
            trans_table = output[5]
            codon_densities = output[6]

            output_line = '%s,%s,%s,%s,%s,%s' % (acc, gc, gc3, trans_table, cds_count, codon_count)
            for frame in codon_densities:
                for stop in sorted(codon_densities[frame]):
                    output_line += ',%s' % (codon_densities[frame][stop])
            output_line += '\n'
            output_file.write(output_line)

    output_file.close()


def write_combined_densities(results):

    output_file = open('outputs/osc_densities/combined_osc_densities_t4.csv', 'w')
    header = 'acc,gc,gc3,trans_table,cds_count,codon_count'
    for frame in frames:
        header += ',%s_density' % (frame)
    header += '\n'
    output_file.write(header)

    for result in results:
        outputs = result.get()
        for output in outputs:
            acc = output[0]
            gc = output[1]
            gc3 = output[2]
            cds_count = output[3]
            codon_count = output[4]
            trans_table = output[5]
            combined_densities = output[7]

            output_line = '%s,%s,%s,%s,%s,%s' % (acc, gc, gc3, trans_table, cds_count, codon_count)
            for frame in combined_densities:
                output_line += ',%s' % (combined_densities[frame])
            output_line += '\n'
            output_file.write(output_line)

    output_file.close()


def write_to_files(results):

    write_codon_densities(results)
    write_combined_densities(results)


def main():


    t0_main = time.time()

    create_directory('outputs/osc_densities/')

    accessions = get_accessions('outputs/genome_extractions_t4/')

    accession_count = 0
    genomes = []
    acc_counts = {}
    for accession in sorted(accessions):
        accession_count += 1
        # if accession == 'AE014184':
        if accession_count:
        # if accession_count <= 1:
            genomes.append(accession)
            acc_counts[accession] = accession_count

    results = run_in_parralell(genomes, [acc_counts], run_genomes)
    # run_genomes(genomes, acc_counts)

    write_to_files(results)





    t1_main = time.time()
    print ('\n%s\nTotal time: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))


if __name__ == '__main__':
    main()
