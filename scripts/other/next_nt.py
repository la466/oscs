#!/usr/bin/python

# Script number:				na
# File:
# Prerequisite script(s):
# Prerequisite file(s):
# Description:				    Get the use of nucleotides that come after an OSC
# Output files:					nt_after_osc



import numpy as np
import re
import time
import multiprocessing as mp
import os
import sys


##########################
# VARIABLES
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

##########################
# FUNCTIONS
##########################


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

def get_stops(cds_seqs):

    stops_counts = {}
    for stop in stop_codons:
        stops_counts[stop] = 0

    for cds in cds_seqs:
        stops_counts[cds[-3:]] += 1

    stop_props = {}
    for stop in stops_counts:
        stop_props[stop] = stops_counts[stop]/sum(stops_counts.values())

    return(stop_props)


def run_genome(accessions, acc_counts, osc_counts):

    outputs = []

    for accession in accessions:

        t0 = time.time()

        print('%s: %s' % (acc_counts[accession], accession))

        cds_seqs = get_seqs(accession)
        stop_props = get_stops(cds_seqs)
        gc, gc3 = calc_gc(cds_seqs)

        osc_props = {}
        for frame in osc_counts[accession]:
            osc_props[frame] = {}
            for stop in osc_counts[accession][frame]:
                osc_props[frame][stop] = osc_counts[accession][frame][stop] / sum(osc_counts[accession][frame].values())

        output = [accession, gc, gc3, stop_props, osc_props]

        outputs.append(output)

        t1 = time.time()
        # print('Genome time: %s' % (t1-t0))

    return (outputs)


def calc_gc(cds_list):

    gc_count = gc1_count = gc2_count = gc3_count = nt_count = nt1_count = nt2_count = nt3_count = 0
    pos_count2 = {}
    for nt in ['A', 'C', 'G', 'T']:
        pos_count2[nt] = 0

    for cds in cds_list:

        nts = list(cds[:-3])
        gc_count += nts.count('G')
        gc_count += nts.count('C')
        nt_count += len(nts)

        nt_index = 0
        for nt in nts:
            nt_index += 1
            if nt_index > 3:
                if nt_index % 3 == 1:
                    nt1_count += 1
                    if nt == 'G' or nt == 'C':
                        gc1_count += 1
                        gc_count += 1
                elif nt_index % 3 == 2:
                    nt2_count += 1
                    pos_count2[nt] += 1
                    if nt == 'G' or nt == 'C':
                        gc2_count += 1
                        gc_count += 1
                else:
                    nt3_count += 1
                    if nt == 'G' or nt == 'C':
                        gc3_count += 1
                        gc_count += 1

    gc = gc_count / float(nt_count)
    gc1 = gc1_count / float(nt1_count)
    gc2 = gc2_count / float(nt2_count)
    gc3 = gc3_count / float(nt3_count)

    return(gc,gc1,gc2,gc3, pos_count2)

def get_accessions(directory):

    accessions = []
    for file in os.listdir(directory):
        if not file.endswith('.DS_Store'):
            accessions.append(re.findall(r'(.+)(?=\.)', file)[0])

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


def get_next_nt(cds_seqs):

    next_nt_1 = {}
    for nt in ['A', 'C', 'G', 'T']:
        next_nt_1[nt] = 0

    for cds in cds_seqs:
        for i in range(1, len(cds)-4, 3):
            if(cds[i:i+3]) in stop_codons:
                next_nt_1[cds[i+4]] += 1

    norm_next_nt_1 = {}
    for nt in next_nt_1:
        norm_next_nt_1[nt] = np.divide(next_nt_1[nt], sum(next_nt_1.values()))

    return(norm_next_nt_1)

def norm_usage_next_nt(norm_next_nt, pos_count2):

    norm_count2 = {}
    for nt in pos_count2:
        norm_count2[nt] = np.divide(pos_count2[nt], sum(pos_count2.values()))

    norm_usage_next = {}
    for nt in norm_next_nt:
        norm_usage_next[nt] = np.divide(norm_next_nt[nt], norm_count2[nt])

    return(norm_usage_next)


def run_genome(accessions, acc_counts):

    outputs = []

    for accession in accessions:
        print('%s: %s' % (acc_counts[accession], accession))

        cds_seqs = get_seqs(accession)
        gc, gc1, gc2, gc3, pos_count2 = calc_gc(cds_seqs)

        norm_next_nt = get_next_nt(cds_seqs)
        norm_usage_next = norm_usage_next_nt(norm_next_nt, pos_count2)
        output = [accession, gc, gc2, norm_usage_next]
        outputs.append(output)

    return (outputs)



def write_to_file(results):

    output_file = open('outputs/other/nt_after_osc.csv', 'w')

    head = 'acc,gc,gc2,a,c,g,t\n'
    output_file.write(head)

    for result in results:
        outputs = result.get()
        for output in outputs:
            acc = output[0]
            gc = output[1]
            gc2 = output[2]
            norm_next_nt = output[3]

            output_line = '%s,%s,%s' % (acc, gc, gc2)
            for nt in sorted(norm_next_nt):
                output_line += ',%s' % (norm_next_nt[nt])
            output_line += '\n'

            output_file.write(output_line)

    output_file.close()


def main():

    t0_main = time.time()

    setupDirectories('outputs/other/')

    accs = get_accessions('outputs/genome_extractions/')

    accessions = []
    acc_counts = {}

    accession_count = 0
    for accession in accs:
        accession_count += 1
        if accession_count:
        # if accession_count <= 1:
            accessions.append(accession)
            acc_counts[accession] = accession_count

    # results = run_genome(accessions, acc_counts)
    results = run_in_parralell(accessions, [acc_counts], run_genome)

    write_to_file(results)

    t1_main = time.time()
    print ('\n%s\nTotal time: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))


if __name__ == '__main__':
    main()
