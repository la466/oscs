#!/usr/bin/python

# Script number:				6.4
# File:							4 of 6
# Prerequisite script(s):
# Prerequisite file(s):
# Description:					Model sequences using fifth order Markov model
# Output files:


import numpy as np
import re
import time
import multiprocessing as mp
import os
import sys
from scipy.stats import t
import scipy.stats
import random

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

stops = {}
stops[4] = ['TAA', 'TAG']
stops[11] = ['TAA', 'TAG', 'TGA']

frames = [1,2,'both']

bases = ['A', 'C', 'G', 'T']

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
            accessions.append(re.findall(r'(.+)(?=\.txt)', file)[0])

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


def get_genome_stops(accession):

    if accession in t4:
        genome_stops = stops[4]
    else:
        genome_stops = stops[11]

    return(genome_stops)


def get_usage(genome_seqs):

    usage = {}
    norm_usage = {}

    for base1 in bases:
        for base2 in bases:
            for base3 in bases:
                for base4 in bases:
                    for base5 in bases:
                        usage[base1+'0'+base2+'1'+base3+'2'+base4+'0'+base5+'1'] = {}
                        usage[base1+'1'+base2+'2'+base3+'0'+base4+'1'+base5+'2'] = {}
                        usage[base1+'2'+base2+'0'+base3+'1'+base4+'2'+base5+'0'] = {}
                        norm_usage[base1+'0'+base2+'1'+base3+'2'+base4+'0'+base5+'1'] = {}
                        norm_usage[base1+'1'+base2+'2'+base3+'0'+base4+'1'+base5+'2'] = {}
                        norm_usage[base1+'2'+base2+'0'+base3+'1'+base4+'2'+base5+'0'] = {}

                        for base6 in bases:
                            usage[base1+'0'+base2+'1'+base3+'2'+base4+'0'+base5+'1'][base6] = 0
                            usage[base1+'1'+base2+'2'+base3+'0'+base4+'1'+base5+'2'][base6] = 0
                            usage[base1+'2'+base2+'0'+base3+'1'+base4+'2'+base5+'0'][base6] = 0
                            norm_usage[base1+'0'+base2+'1'+base3+'2'+base4+'0'+base5+'1'][base6] = 0
                            norm_usage[base1+'1'+base2+'2'+base3+'0'+base4+'1'+base5+'2'][base6] = 0
                            norm_usage[base1+'2'+base2+'0'+base3+'1'+base4+'2'+base5+'0'][base6] = 0

    count = 0
    for cds in genome_seqs:
        count += 1
        if count:

            query_seq = cds[3:-3]

            # Use the first 5 nts of the query as starter seed sequence
            for i in range(5,len(query_seq)-5):
                pos1 = query_seq[i]
                pos2 = query_seq[i+1]
                pos3 = query_seq[i+2]
                pos4 = query_seq[i+3]
                pos5 = query_seq[i+4]
                pos6 = query_seq[i+5]

                index1 = str(i%3)
                index2 = str((i+1)%3)
                index3 = str((i+2)%3)
                index4 = str(i%3)
                index5 = str((i+1)%3)

                usage[pos1+index1+pos2+index2+pos3+index3+pos4+index4+pos5+index5][pos6] += 1


    for pair in sorted(usage):
        for base in sorted(usage[pair]):
            norm_usage[pair][base] = np.divide(usage[pair][base],sum(usage[pair].values()))

    return(norm_usage)


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


def get_off_frame_densities(seqs):

    off_frame_frequencies = {}

    for frame in frames:
        off_frame_frequencies[frame] = {}
        for codon in codon_list:
            off_frame_frequencies[frame][codon] = 0

    codon_count = 0
    for cds in seqs:

        codon_count += len(cds[6:-3])/3

        f1 = re.findall(r'.{3}', cds[7:])
        f2 = re.findall(r'.{3}', cds[8:])
        for codon in codon_list:
            off_frame_frequencies[1][codon] += f1.count(codon)
            off_frame_frequencies[2][codon] += f2.count(codon)

    off_frame_densities = {}
    for frame in off_frame_frequencies:
        off_frame_densities[frame] = {}
        for codon in off_frame_frequencies[frame]:
            off_frame_densities[frame][codon] = (100/codon_count)*off_frame_frequencies[frame][codon]

    return(off_frame_densities, codon_count)




def generate_seqs(repeat_list, genome_seqs, norm_usage):


    samples = []

    for repeat in repeat_list:

        np.random.seed()

        print('Repeat %s' % (repeat+1))

        generated_seqs = []

        gen_count = 0
        for cds in genome_seqs:
            gen_count += 1
            if gen_count:

                new_seq = cds[0:8]
                required_nts = len(cds) - 8 - 3  # minus start seq and stop codon

                pos1 = cds[3]
                pos2 = cds[4]
                pos3 = cds[5]
                pos4 = cds[6]
                pos5 = cds[7]

                index1 = '0'
                index2 = '1'
                index3 = '2'
                index4 = '0'
                index5 = '1'

                for i in range(1,required_nts+1):
                    cum_prob = 0

                    prob = random.random()

                    pair = pos1+index1+pos2+index2+pos3+index3+pos4+index4+pos5+index5
                    for base in norm_usage[pair]:
                        cum_prob += norm_usage[pair][base]

                        if cum_prob > prob:
                            new_seq += base

                            pos1 = pos2
                            pos2 = pos3
                            pos3 = pos4
                            pos4 = pos5
                            pos5 = base


                            index1 = index2
                            index2 = index3
                            index3 = index4
                            index4 = index5
                            index5 = index2

                            break

                new_seq += cds[-3:]


                generated_seqs.append(new_seq)

        gen_off_frame_densities, codon_count = get_off_frame_densities(generated_seqs)
        gc, gc3 = calc_gc(generated_seqs)

        sample = [repeat+1, gc, gc3, codon_count, gen_off_frame_densities]
        samples.append(sample)


    return(samples)

def write_to_file(accession, real_sample, results):

    output_file = open('outputs/mmodels/53mm/' + accession + '_53mm.csv', 'w')
    header_string = 'rep,gc,gc3,codon_count'
    for frame in [1,2]:
        for codon in sorted(codon_list):
            header_string += ',%s_%s' % (codon, frame)
    header_string += '\n'
    output_file.write(header_string)

    write_line(accession, real_sample, output_file)

    for result in results:
        outputs = result.get()
        for output in outputs:
            write_line(accession, output, output_file)

    output_file.close()


def write_line(accession, sample, output_file):

    repeat = sample[0]

    gc = sample[1]
    gc3 = sample[2]
    codon_count = sample[3]
    codon_densities = sample[4]

    output_line = '%s_%s,%s,%s,%s' % (accession, repeat, gc, gc3, codon_count)
    for frame in [1,2]:
        for codon in sorted(codon_densities[frame]):
            output_line += ',%s' % (codon_densities[frame][codon])
    output_line += '\n'
    output_file.write(output_line)




def run_genomes(genomes, acc_counts):

    required_repeats = 200

    for accession in genomes:

        print('%s: %s' % (acc_counts[accession], accession))

        genome_stops = get_genome_stops(accession)
        genome_seqs = get_seqs(accession)
        gc, gc3 = calc_gc(genome_seqs)
        real_densities, real_codon_count = get_off_frame_densities(genome_seqs)
        real_sample = ['real', gc, gc3, real_codon_count, real_densities]

        norm_usage = get_usage(genome_seqs)

        repeat_list = list(range(required_repeats))
        # results = generate_seqs(repeat_list, genome_seqs, norm_usage)
        results = run_in_parralell(repeat_list, [genome_seqs, norm_usage], generate_seqs)

        write_to_file(accession, real_sample, results)




def main():

    t0_main = time.time()

    create_directory('outputs/mmodels/')
    create_directory('outputs/mmodels/53mm/')

    accessions = get_accessions('outputs/genome_extractions/')

    genomes = []
    acc_counts = {}

    accession_count = 0
    for accession in sorted(accessions):
        accession_count += 1
        if accession_count:
        # if accession_count <= 1:
            genomes.append(accession)
            acc_counts[accession] = accession_count


    run_genomes(genomes, acc_counts)
    # results = run_in_parralell(genomes, [acc_counts], run_genome)





    t1_main = time.time()
    print ('\n%s\nTotal time: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))


if __name__ == '__main__':
    main()
