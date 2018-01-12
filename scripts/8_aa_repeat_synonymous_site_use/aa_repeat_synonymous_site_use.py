#!/usr/bin/python

# Script number:				8.1
# File:							1 of 1
# Prerequisite script(s):
# Prerequisite file(s):
# Description:					Get the use of nucleotides at synonymous sites for amino acids whos codons when repeated can encode an OSC, followed by another codon which strictly prevents an OSC
# Output files:                 aa_repeat_osc_synonymous_sites.csv, aa_repeat_osc_synonymous_sites_only_at.csv


import numpy as np
import re
import time
import multiprocessing as mp
import os
import sys
import random
import collections


##########################
# VARIABLES
##########################

codon_map = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"St", "TCC":"St", "TCA":"St", "TCG":"St",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
    "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"Sa", "AGC":"Sa", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}


codon_list = [codon for codon in sorted(codon_map)]

frames = [1,-1]
nts = ['A', 'C', 'G', 'T']


##########################
# FUNCTIONS #
##########################


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






def analyse(cds_seqs):

    cds_count = 0

    nt_use = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict((int))))

    nt_use_at_gc = {}
    for aa in ['I', 'V']:
        nt_use_at_gc[aa] = {}
        for site in [3,6]:
            nt_use_at_gc[aa][site] = {}
            for nt in ['A', 'T']:
                nt_use_at_gc[aa][site][nt] = 0

    pair_count = 0
    i_at_count = 0

    for cds in cds_seqs:
        cds_count += 1
        if cds_count:

            # Get a list of codons and amino acids
            codons = re.findall(r'.{3}', cds[3:-3])
            aas = [codon_map[codon] for codon in codons]

            for i in range(len(aas)-2):

                pair_count += 1

                # For +1 TAA
                if aas[i] == 'I' and aas[i+1] == 'I' and codons[i+2][0] in ['C', 'T']:
                    nt_use[aas[i]][3][codons[i][-1]] += 1
                    nt_use[aas[i]][6][codons[i+1][-1]] += 1

                    if codons[i][-1] in ['A', 'T'] and codons[i+1][-1] in ['A', 'T']:

                        i_at_count += 1

                        if codons[i][-1] == 'A':
                            nt_use_at_gc['I'][3]['A'] += 1
                        else:
                            nt_use_at_gc['I'][3]['T'] += 1
                        if codons[i+1][-1] == 'A':
                            nt_use_at_gc['I'][6]['A'] += 1
                        else:
                            nt_use_at_gc['I'][6]['T'] += 1

                # For +1 TAG
                if aas[i] == 'V' and aas[i+1] == 'V' and codons[i+2][0] in ['C', 'T']:
                    nt_use[aas[i]][3][codons[i][-1]] += 1
                    nt_use[aas[i]][6][codons[i+1][-1]] += 1

                    if codons[i][-1] in ['A', 'T'] and codons[i+1][-1] in ['A', 'T']:
                        if codons[i][-1] == 'A':
                            nt_use_at_gc['V'][3]['A'] += 1
                        else:
                            nt_use_at_gc['V'][3]['T'] += 1
                        if codons[i+1][-1] == 'A':
                            nt_use_at_gc['V'][6]['A'] += 1
                        else:
                            nt_use_at_gc['V'][6]['T'] += 1


                # For +2 TAA
                if aas[i] == 'N' and aas[i+1] == 'N' and codons[i+2][:2] not in ['AA', 'AG', 'GA']:
                    nt_use[aas[i]][3][codons[i][-1]] += 1
                    nt_use[aas[i]][6][codons[i+1][-1]] += 1

                # For +2 TAG
                if aas[i] == 'Sa' and aas[i+1] == 'Sa' and codons[i+2][:2] not in ['AA', 'AG', 'GA']:
                    nt_use[aas[i]][3][codons[i][-1]] += 1
                    nt_use[aas[i]][6][codons[i+1][-1]] += 1

                # For +2 TGA
                if aas[i] == 'D' and aas[i+1] == 'D' and codons[i+2][:2] not in ['AA', 'AG', 'GA']:
                    nt_use[aas[i]][3][codons[i][-1]] += 1
                    nt_use[aas[i]][6][codons[i+1][-1]] += 1


    norm_nt_use = {}
    for aa in sorted(nt_use):
        norm_nt_use[aa] = {}
        for site in sorted(nt_use[aa]):
            norm_nt_use[aa][site] = {}
            for nt in nts:
                if nt in nt_use[aa][site]:
                    norm_nt_use[aa][site][nt] = nt_use[aa][site][nt] / sum(nt_use[aa][site].values())
                else:
                    norm_nt_use[aa][site][nt] = 0

    norm_nt_use_at_gc = {}
    for aa in sorted(nt_use_at_gc):
        norm_nt_use_at_gc[aa] = {}
        for site in sorted(nt_use_at_gc[aa]):
            norm_nt_use_at_gc[aa][site] = {}
            for nt in sorted(nt_use_at_gc[aa][site]):
                if sum(nt_use_at_gc[aa][site].values()) == 0:
                    norm_nt_use_at_gc[aa][site][nt] = 0
                else:
                    norm_nt_use_at_gc[aa][site][nt] = nt_use_at_gc[aa][site][nt] / sum(nt_use_at_gc[aa][site].values())



    return(norm_nt_use, norm_nt_use_at_gc, pair_count, i_at_count)


def run_genome(accessions, acc_counts):

    outputs = []

    for accession in accessions:


        print('Genome %s: %s' % (acc_counts[accession], accession))

        cds_seqs = get_seqs(accession)
        gc, gc3 = calc_gc(cds_seqs)

        norm_nt_use, norm_nt_use_at_gc, pair_count, i_at_count = analyse(cds_seqs)
        output = [accession, gc, gc3, norm_nt_use, norm_nt_use_at_gc, pair_count, i_at_count]
        outputs.append(output)

    return(outputs)


def write_results_all(results):

    output_file = open('outputs/aa_repeat_synonymous_sites/aa_repeat_osc_synonymous_sites.csv', 'w')
    output_file.write('acc,gc,gc3,TAA_1_a3,TAA_1_a6,TAG_1_a3,TAG_1_a6,TAA_2_t3,TAA_2_t6,TAG_2_t3,TAG_2_t6,TGA_2_t3,TGA_2_t6\n')

    for result in results:
        outputs = result.get()
        for output in outputs:
            acc = output[0]
            gc = output[1]
            gc3 = output[2]
            norm_nt_use = output[3]

            output_line = '%s,%s,%s' % (acc,gc,gc3)
            output_line += ',%s,%s' % (norm_nt_use['I'][3]['A'], norm_nt_use['I'][6]['A'])
            output_line += ',%s,%s' % (norm_nt_use['V'][3]['A'], norm_nt_use['V'][6]['A'])
            output_line += ',%s,%s' % (norm_nt_use['N'][3]['T'], norm_nt_use['N'][6]['T'])
            output_line += ',%s,%s' % (norm_nt_use['Sa'][3]['T'], norm_nt_use['Sa'][6]['T'])
            output_line += ',%s,%s' % (norm_nt_use['D'][3]['T'], norm_nt_use['D'][6]['T'])
            output_line += '\n'

            output_file.write(output_line)

    output_file.close()

def write_results_at_gc(results):

    output_file = open('outputs/aa_repeat_synonymous_sites/aa_repeat_osc_synonymous_sites_only_at.csv', 'w')
    output_file.write('acc,gc,gc3,TAA_1_a3,TAA_1_a6,TAA_1_t3,TAA_1_t6,TAG_1_a3,TAG_1_a6,TAG_1_t3,TAG_1_t6\n')

    for result in results:
        outputs = result.get()
        for output in outputs:
            acc = output[0]
            gc = output[1]
            gc3 = output[2]
            norm_nt_use_at_gc = output[4]

            output_line = '%s,%s,%s' % (acc,gc,gc3)
            output_line += ',%s,%s' % (norm_nt_use_at_gc['I'][3]['A'], norm_nt_use_at_gc['I'][6]['A'])
            output_line += ',%s,%s' % (norm_nt_use_at_gc['I'][3]['T'], norm_nt_use_at_gc['I'][6]['T'])
            output_line += ',%s,%s' % (norm_nt_use_at_gc['V'][3]['A'], norm_nt_use_at_gc['V'][6]['A'])
            output_line += ',%s,%s' % (norm_nt_use_at_gc['V'][3]['T'], norm_nt_use_at_gc['V'][6]['T'])
            output_line += '\n'

            output_file.write(output_line)

    output_file.close()


def main():


    t0_main = time.time()

    create_directory('outputs/aa_repeat_synonymous_sites/')
    accessions = get_accessions('outputs/genome_extractions/')

    genomes = []
    acc_counts = {}

    accession_count = 0
    for accession in sorted(accessions):
        accession_count += 1
        # if accession == 'CP003881':
        # if accession == 'CP001631':
        # if accession_count >= 77 and accession_count < 80:
        # if accession_count <= 10:
        if accession_count:
            genomes.append(accession)
            acc_counts[accession] = accession_count

    # run_genome(genomes, acc_counts)
    results = run_in_parralell(genomes, [acc_counts], run_genome)
    write_results_all(results)
    write_results_at_gc(results)

    total_pairs = 0
    total_i_at = 0
    for result in results:
        outputs = result.get()
        for output in outputs:
            pair_count = output[5]
            i_at_count = output[6]

            total_pairs += pair_count
            total_i_at += i_at_count

    print(total_i_at)
    print(total_pairs)
    print(total_i_at/total_pairs)


    t1_main = time.time()
    print ('\n%s\nTotal time: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))


if __name__ == '__main__':
    main()
