#!/usr/bin/python

# Script number:			    7.
# File:
# Prerequisite script(s):
# Prerequisite file(s):
# Description:					Randomise synonymous sites within the genome, preserving amino acid identities and genome codon usage biases. Outputs OSC densities for 200 simulations for each genome. Codon count excludes one-fold degenerates for t4 genomes
# Output files:


import numpy as np
import re
import time
import multiprocessing as mp
import os
import sys
import scipy.stats
import random

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--acc", help="Select an accession to run")
parser.add_argument("--start", help="Select a genome number to start the run")
parser.add_argument("--end", help="Select a genome number to end the run")
args = parser.parse_args()

##########################
# Variables
##########################


codon_map = {}
codon_map[11] = {"TTT":"F", "TTC":"F", "TTA":"Lt", "TTG":"Lt",
    "TCT":"St", "TCC":"St", "TCA":"St", "TCG":"St",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
    "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
    "CTT":"Lc", "CTC":"Lc", "CTA":"Lc", "CTG":"Lc",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"Rc", "CGC":"Rc", "CGA":"Rc", "CGG":"Rc",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"Sa", "AGC":"Sa", "AGA":"Ra", "AGG":"Ra",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

codon_map[4] = {"TTT":"F", "TTC":"F", "TTA":"Lt", "TTG":"Lt",
    "TCT":"St", "TCC":"St", "TCA":"St", "TCG":"St",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
    "TGT":"C", "TGC":"C", "TGA":"W", "TGG":"W",
    "CTT":"Lc", "CTC":"Lc", "CTA":"Lc", "CTG":"Lc",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"Rc", "CGC":"Rc", "CGA":"Rc", "CGG":"Rc",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"Sa", "AGC":"Sa", "AGA":"Ra", "AGG":"Ra",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}


synonymous_site_blocks = {}
synonymous_site_blocks[11] = {
    'TTT': ['TTT', 'TTC'],
    'TTC': ['TTT', 'TTC'],
    'TTA': ['TTA', 'TTG'],
    'TTG': ['TTA', 'TTG'],
    'CTT': ['CTT', 'CTC', 'CTA', 'CTG'],
    'CTC': ['CTT', 'CTC', 'CTA', 'CTG'],
    'CTA': ['CTT', 'CTC', 'CTA', 'CTG'],
    'CTG': ['CTT', 'CTC', 'CTA', 'CTG'],
    'ATT': ['ATT', 'ATC', 'ATA'],
    'ATC': ['ATT', 'ATC', 'ATA'],
    'ATA': ['ATT', 'ATC', 'ATA'],
    'ATG': ['ATG'],
    'GTT': ['GTT', 'GTC', 'GTA', 'GTG'],
    'GTC': ['GTT', 'GTC', 'GTA', 'GTG'],
    'GTA': ['GTT', 'GTC', 'GTA', 'GTG'],
    'GTG': ['GTT', 'GTC', 'GTA', 'GTG'],
    'TCT': ['TCT', 'TCC', 'TCA', 'TCG'],
    'TCC': ['TCT', 'TCC', 'TCA', 'TCG'],
    'TCA': ['TCT', 'TCC', 'TCA', 'TCG'],
    'TCG': ['TCT', 'TCC', 'TCA', 'TCG'],
    'CCT': ['CCT', 'CCC', 'CCA', 'CCG'],
    'CCC': ['CCT', 'CCC', 'CCA', 'CCG'],
    'CCA': ['CCT', 'CCC', 'CCA', 'CCG'],
    'CCG': ['CCT', 'CCC', 'CCA', 'CCG'],
    'ACT': ['ACT', 'ACC', 'ACA', 'ACG'],
    'ACC': ['ACT', 'ACC', 'ACA', 'ACG'],
    'ACA': ['ACT', 'ACC', 'ACA', 'ACG'],
    'ACG': ['ACT', 'ACC', 'ACA', 'ACG'],
    'GCT': ['GCT', 'GCC', 'GCA', 'GCG'],
    'GCC': ['GCT', 'GCC', 'GCA', 'GCG'],
    'GCA': ['GCT', 'GCC', 'GCA', 'GCG'],
    'GCG': ['GCT', 'GCC', 'GCA', 'GCG'],
    'TAT': ['TAT', 'TAC'],
    'TAC': ['TAT', 'TAC'],
    'CAT': ['CAT', 'CAC'],
    'CAC': ['CAT', 'CAC'],
    'CAA': ['CAA', 'CAG'],
    'CAG': ['CAA', 'CAG'],
    'AAT': ['AAT', 'AAC'],
    'AAC': ['AAT', 'AAC'],
    'AAA': ['AAA', 'AAG'],
    'AAG': ['AAA', 'AAG'],
    'GAT': ['GAT', 'GAC'],
    'GAC': ['GAT', 'GAC'],
    'GAA': ['GAA', 'GAG'],
    'GAG': ['GAA', 'GAG'],
    'TGT': ['TGT', 'TGC'],
    'TGC': ['TGT', 'TGC'],
    'TGG': ['TGG'],
    'CGT': ['CGT', 'CGC', 'CGA', 'CGG'],
    'CGC': ['CGT', 'CGC', 'CGA', 'CGG'],
    'CGA': ['CGT', 'CGC', 'CGA', 'CGG'],
    'CGG': ['CGT', 'CGC', 'CGA', 'CGG'],
    'AGT': ['AGT', 'AGC'],
    'AGC': ['AGT', 'AGC'],
    'AGA': ['AGA', 'AGG'],
    'AGG': ['AGA', 'AGG'],
    'GGT': ['GGT', 'GGC', 'GGA', 'GGG'],
    'GGC': ['GGT', 'GGC', 'GGA', 'GGG'],
    'GGA': ['GGT', 'GGC', 'GGA', 'GGG'],
    'GGG': ['GGT', 'GGC', 'GGA', 'GGG'],
    'TAA': [],
    'TAG': [],
    'TGA': []
}
synonymous_site_blocks[4] = {
    'TTT': ['TTT', 'TTC'],
    'TTC': ['TTT', 'TTC'],
    'TTA': ['TTA', 'TTG'],
    'TTG': ['TTA', 'TTG'],
    'CTT': ['CTT', 'CTC', 'CTA', 'CTG'],
    'CTC': ['CTT', 'CTC', 'CTA', 'CTG'],
    'CTA': ['CTT', 'CTC', 'CTA', 'CTG'],
    'CTG': ['CTT', 'CTC', 'CTA', 'CTG'],
    'ATT': ['ATT', 'ATC', 'ATA'],
    'ATC': ['ATT', 'ATC', 'ATA'],
    'ATA': ['ATT', 'ATC', 'ATA'],
    'ATG': ['ATG'],
    'GTT': ['GTT', 'GTC', 'GTA', 'GTG'],
    'GTC': ['GTT', 'GTC', 'GTA', 'GTG'],
    'GTA': ['GTT', 'GTC', 'GTA', 'GTG'],
    'GTG': ['GTT', 'GTC', 'GTA', 'GTG'],
    'TCT': ['TCT', 'TCC', 'TCA', 'TCG'],
    'TCC': ['TCT', 'TCC', 'TCA', 'TCG'],
    'TCA': ['TCT', 'TCC', 'TCA', 'TCG'],
    'TCG': ['TCT', 'TCC', 'TCA', 'TCG'],
    'CCT': ['CCT', 'CCC', 'CCA', 'CCG'],
    'CCC': ['CCT', 'CCC', 'CCA', 'CCG'],
    'CCA': ['CCT', 'CCC', 'CCA', 'CCG'],
    'CCG': ['CCT', 'CCC', 'CCA', 'CCG'],
    'ACT': ['ACT', 'ACC', 'ACA', 'ACG'],
    'ACC': ['ACT', 'ACC', 'ACA', 'ACG'],
    'ACA': ['ACT', 'ACC', 'ACA', 'ACG'],
    'ACG': ['ACT', 'ACC', 'ACA', 'ACG'],
    'GCT': ['GCT', 'GCC', 'GCA', 'GCG'],
    'GCC': ['GCT', 'GCC', 'GCA', 'GCG'],
    'GCA': ['GCT', 'GCC', 'GCA', 'GCG'],
    'GCG': ['GCT', 'GCC', 'GCA', 'GCG'],
    'TAT': ['TAT', 'TAC'],
    'TAC': ['TAT', 'TAC'],
    'CAT': ['CAT', 'CAC'],
    'CAC': ['CAT', 'CAC'],
    'CAA': ['CAA', 'CAG'],
    'CAG': ['CAA', 'CAG'],
    'AAT': ['AAT', 'AAC'],
    'AAC': ['AAT', 'AAC'],
    'AAA': ['AAA', 'AAG'],
    'AAG': ['AAA', 'AAG'],
    'GAT': ['GAT', 'GAC'],
    'GAC': ['GAT', 'GAC'],
    'GAA': ['GAA', 'GAG'],
    'GAG': ['GAA', 'GAG'],
    'TGT': ['TGT', 'TGC'],
    'TGC': ['TGT', 'TGC'],
    'TGG': ['TGG', 'TGA'],
    'CGT': ['CGT', 'CGC', 'CGA', 'CGG'],
    'CGC': ['CGT', 'CGC', 'CGA', 'CGG'],
    'CGA': ['CGT', 'CGC', 'CGA', 'CGG'],
    'CGG': ['CGT', 'CGC', 'CGA', 'CGG'],
    'AGT': ['AGT', 'AGC'],
    'AGC': ['AGT', 'AGC'],
    'AGA': ['AGA', 'AGG'],
    'AGG': ['AGA', 'AGG'],
    'GGT': ['GGT', 'GGC', 'GGA', 'GGG'],
    'GGC': ['GGT', 'GGC', 'GGA', 'GGG'],
    'GGA': ['GGT', 'GGC', 'GGA', 'GGG'],
    'GGG': ['GGT', 'GGC', 'GGA', 'GGG'],
    'TAA': [],
    'TAG': [],
    'TGA': ['TGA', 'TGG']
}


codon_list = [codon for codon in sorted(codon_map[11])]

bases = ['A', 'C', 'G', 'T']

stops = {}
stops[4]= ['TAA', 'TAG']
stops[11] = ['TAA', 'TAG', 'TGA']

frames = [1, 2]

required_repeats = 200

# def get_t4_genomes():
#
# 	t4_genomes = []
# 	with open('outputs/gene_filtering/table_4_genomes_in_list.txt', 'rU') as t4_file:
# 		lines = t4_file.readlines()
# 		for line in lines:
# 			acc = line.strip('\n')
# 			t4_genomes.append(acc)
# 	return(t4_genomes)
#
# t4_genomes = get_t4_genomes()




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



def get_off_frame_frequencies(accession, seqs):

# if accession in t4_genomes:
    onefold = ['ATG']
    # else:
    #     onefold = ['ATG', 'TGG']


    off_frame_frequencies = {}
    codon_count = 0


    for frame in frames:
        off_frame_frequencies[frame] = {}
        for codon in codon_list:
            off_frame_frequencies[frame][codon] = 0

    for cds in seqs:

        for i in range(0, len(cds)-3, 3):

            codon = cds[i:i+3]

            if codon not in onefold:
                codon_count += 1
                off_frame_frequencies[1][cds[i+1:i+4]] += 1
                off_frame_frequencies[2][cds[i+2:i+5]] += 1

    return(off_frame_frequencies, codon_count)


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


def genome_stats(accession, cds_seqs):

    stats = []
    gc, gc3 = calc_gc(cds_seqs)
    off_frame_counts, codon_count = get_off_frame_frequencies(accession, cds_seqs)

    stats = [gc, gc3, off_frame_counts, codon_count]

    return(stats)



# Write the file header line
def write_header(output_file):

    line = 'acc,gc,gc3,codon_count'

    for frame in sorted(frames):
        for codon in sorted(codon_list):
            line += ',%s_%s' % (codon, frame)

    line += '\n'
    output_file.write(line)


# Write the real results to file
def write_real(output_file, accession, stats):

    gc = stats[0]
    gc3 = stats[1]
    off_frame_counts = stats[2]
    codon_count = stats[3]


    line = '%s,%s,%s,%s' % (accession, gc, gc3,codon_count)

    for frame in sorted(off_frame_counts):
        for codon in sorted(off_frame_counts[frame]):
            line += ',%s' % ((100/codon_count)*off_frame_counts[frame][codon])

    line += '\n'
    output_file.write(line)

# Write the repeat results to file
def write_repeats(output_file, accession, repeat, stats):

    # For each of the repeats, write the line
    # for repeat in stats:

    gc = stats[0]
    gc3 = stats[1]
    off_frame_counts = stats[2]
    codon_count = stats[3]


    line = '%s_r%s,%s,%s,%s' % (accession, repeat+1, gc, gc3, codon_count)

    for frame in sorted(off_frame_counts):
        for codon in sorted(off_frame_counts[frame]):
            line += ',%s' % ((100/codon_count)*off_frame_counts[frame][codon])

    line += '\n'
    output_file.write(line)



def get_synonymous_sites(cds_seqs, table):

    synonymous_sites = {}
    for codon in codon_map[table]:
        if codon_map[table][codon] not in synonymous_sites:
            synonymous_sites[codon_map[table][codon]] = {}
            for base in bases:
                synonymous_sites[codon_map[table][codon]][base] = 0


    for cds in cds_seqs:
        for i in range(3, len(cds)-3, 3):
            # synonymous_usage[cds[i]] += 1
            synonymous_sites[codon_map[table][cds[i:i+3]]][cds[i+2]] += 1

    norm_sites = {}
    for aa in synonymous_sites:
        norm_sites[aa] = {}
        for base in synonymous_sites[aa]:
            norm_sites[aa][base] = np.divide(synonymous_sites[aa][base], sum(synonymous_sites[aa].values()))

    return(norm_sites)





def get_cds_codons(cds_seqs):

    cds_codons = []
    for cds in cds_seqs:
        codons = re.findall(r'.{3}', cds)
        cds_codons.append(codons)

    return(cds_codons)


def randomise(repeats, accession, cds_seqs, synonymous_sites, genome_stops, table):

    outputs = []

    for repeat in repeats:

        print('Repeat %s' % (repeat+1))

        np.random.seed()

        new_seqs = []

        cds_count = 0
        for cds in cds_seqs:
            cds_count += 1
            new_seq = cds[:3]

            for i in range(3,len(cds)-3,3):
                cum_prob = 0

                prob = random.random()
                codon_start = cds[i:i+2]
                new_seq += codon_start

                for base in synonymous_sites[codon_map[table][cds[i:i+3]]]:
                    cum_prob += synonymous_sites[codon_map[table][cds[i:i+3]]][base]

                    if cum_prob > prob:
                        new_seq += base
                        break

                # if str(cds[i:i+2]) + str(random.choice(synonymous_sites)) not in genome_stops:
                #     new_seq += str(cds[i:i+2]) + str(random.choice(synonymous_sites))


            new_seq += cds[-3:]


            new_seqs.append(new_seq)

        stats = genome_stats(accession, new_seqs)

        # print(stats)

        output = [accession, repeat, stats]
        outputs.append(output)


    return(outputs)


def get_table(accession):

    # table = 11
    # if accession in t4_genomes:
    table = 4

    return(table)


def run_genome(accession, acc_count):


    print('%s: %s' % (acc_count, accession))

    cds_seqs = get_seqs(accession)
    table = get_table(accession)

    # Get the usage of each codon compared with synonyms
    synonymous_sites = get_synonymous_sites(cds_seqs, table)

    # Get the codon lists
    cds_codons = get_cds_codons(cds_seqs)


    real_stats = genome_stats(accession, cds_seqs)


    # if accession in t4_genomes:
    genome_stops = stops[4]
    # else:
    #     genome_stops = stops[11]


    results = []
    repeats = list(range(required_repeats))
    # results = randomise(repeats, accession, cds_seqs, synonymous_sites)

    results = run_in_parralell(repeats, [accession, cds_seqs, synonymous_sites, genome_stops, table], randomise)

    output_file = open('outputs/simulation_synonymous_site_t4/%s_synonymous_site.csv' % (accession), 'w')
    write_header(output_file)
    write_real(output_file, accession, real_stats)

    # for output in results:



    for result in results:
        outputs = result.get()
        for output in outputs:
            accession = output[0]
            repeat = output[1]
            repeat_stats = output[2]

            write_repeats(output_file, accession, repeat, repeat_stats)


    output_file.close()


def main_run(accession, accession_count):

    t0 = time.time()

    run_genome(accession, accession_count)

    t1 = time.time()
    print('\nGenome time: %s' % (t1-t0))
    print('-' * 30)



def main():


    t0_main = time.time()

    # Setup output directory
    create_directory('outputs/simulation_synonymous_site_t4/')

    # Get a list of accessions
    accessions = get_accessions('outputs/genome_extractions_t4/')

    print('\nRandomising synonymous sites in CDS for %s genomes\n' % len(accessions))



    # Misc setup
    genomes = []
    acc_counts = {}
    accession_count = 0

    for accession in sorted(accessions):
        accession_count += 1

        if args.acc:
            if accession == args.acc:
                main_run(accession, accession_count)

        elif args.start:
            start = int(args.start)
            if args.end:
                end = int(args.end)
                if accession_count >= start and accession_count < end:
                    main_run(accession, accession_count)
            else:
                if accession_count >= start:
                    main_run(accession, accession_count)

        elif args.end:
            end = int(args.end)
            if accession_count < end:
                main_run(accession, accession_count)
        else:
            main_run(accession, accession_count)



    t1_main = time.time()
    print ('\n%s\nTotal time: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))


if __name__ == '__main__':
    main()
