#!/usr/bin/python

# Script number:				5.1
# File:							1 of 1
# Prerequisite script(s):
# Prerequisite file(s):
# Description:					Repeat the analysis by Seligmann and Pollock 2004 on genome sample.
# Output files:					correlation.csv

import numpy as np
import re
import time
import multiprocessing as mp
import os
import sys
import scipy.stats


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

stops = ['TAA', 'TAG', 'TGA']

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



def calc_gc(cds_list):

    gc_count = nt_count = 0
    for cds in cds_list:

        nts = list(cds[:-3])
        gc_count += nts.count('G')
        gc_count += nts.count('C')
        nt_count += len(nts)

    return gc_count / float(nt_count)


def codon_contribution(codon_list):

    codon_contributions = {}

    for codon in codon_list:

        if codon not in stops:

            contribution = 0

            # Get contribution to +1 fs
            for base in bases:
                fs1 = codon[1:] + base
                if fs1 in stops:
                    contribution += 1

            # Get contribution to +2 fs
            for base1 in bases:
                for base2 in bases:
                    fs2 = codon[-1]+base1+base2
                    if fs2 in stops:
                        contribution += 1

            # Get contribution to -1 fs
            for base in bases:
                fsm1 = base+codon[:2]
                if fsm1 in stops:
                    contribution += 1

            # Get contribution to -2 fs
            for base1 in bases:
                for base2 in bases:
                    fsm2 = base1+base2+codon[0]
                    if fsm2 in stops:
                        contribution += 1

            if contribution not in codon_contributions:
                codon_contributions[contribution] = [codon]
            else:
                codon_contributions[contribution].append(codon)

    return(codon_contributions)


def get_codon_use(cds_seqs):

    codon_use = {}
    for codon in codon_list:
        codon_use[codon] = 0

    codon_count = 0

    for cds in cds_seqs:

        codons = re.findall(r'.{3}', cds[3:-3])

        codon_count += len(codons)

        for codon in codon_use:
            codon_use[codon] += codons.count(codon)


    return(codon_use, codon_count)



def compare_use(codon_use, codon_count, codon_contributions):


    contributions = []
    mean_usages = []
    total_usages = []


    for contribution in sorted(codon_contributions):

        codons = codon_contributions[contribution]

        contribution_codon_use_count = 0
        for codon in codons:
            contribution_codon_use_count += codon_use[codon]
        total_usage = contribution_codon_use_count / codon_count

        perc_usages = []
        for codon in codons:
            codon_usage = codon_use[codon]
            perc_usages.append(codon_usage/codon_count)

        mean_usage = np.mean(perc_usages)

        contributions.append(contribution)
        mean_usages.append(mean_usage)
        total_usages.append(total_usage)


    corr_mean = scipy.stats.pearsonr(contributions, mean_usages)
    corr_total = scipy.stats.pearsonr(contributions, total_usages)

    return(corr_mean, corr_total)



def run_genome(accessions, acc_counts, codon_contributions):

    outputs = []

    for accession in accessions:


        print('%s: %s' % (acc_counts[accession], accession))

        cds_seqs = get_seqs(accession)

        gc = calc_gc(cds_seqs)

        codon_use, codon_count = get_codon_use(cds_seqs)

        corr_mean, corr_total = compare_use(codon_use, codon_count, codon_contributions)


        output = [accession, gc, corr_mean, corr_total]
        outputs.append(output)


    return(outputs)



def main():


    t0_main = time.time()

    accessions = get_accessions('outputs/genome_extractions/')
    create_directory('outputs/seligmann_pollock_repeat/')

    codon_contributions = codon_contribution(codon_list)

    genomes = []
    acc_counts = {}

    accession_count = 0
    for accession in accessions:
        accession_count += 1
        if accession_count:
        # if accession_count < 10:
        # if accession == 'CP000013':
            genomes.append(accession)
            acc_counts[accession] = accession_count

    results = run_in_parralell(genomes, [acc_counts, codon_contributions], run_genome)


    output_file = open('outputs/seligmann_pollock_repeat/correlation.csv', 'w')
    output_file.write('acc,gc,cor_mean_use,p_mean_use,cor_total_use,cor_total_p\n')

    for result in results:
        outputs = result.get()
        for output in outputs:
            acc = output[0]
            gc = output[1]
            corr_mean = output[2]
            corr_total = output[3]

            output_line = '%s,%s,%s,%s,%s,%s\n' % (acc, gc, corr_mean[0], corr_mean[1],corr_total[0],corr_total[1])
            output_file.write(output_line)

    output_file.close()




    t1_main = time.time()
    print ('\n%s\nTotal time: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))


if __name__ == '__main__':
    main()
