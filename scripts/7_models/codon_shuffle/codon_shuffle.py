#!/usr/bin/python

# Script number:			    7.1
# File:							1 of 9
# Prerequisite script(s):
# Prerequisite file(s):
# Description:					Randomise CDSs by shuffling the codons within each CDS. Outputs OSC densities for 200 simulations for each genome.
# Output files:



import numpy as np
import re
import time
import multiprocessing as mp
import os
import sys
import itertools

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
frames = [1, 2]
required_repeats = 200


##########################
# FUNCTIONS
##########################


def create_directory(directory_path):

	if os.path.exists(directory_path):
		print ('Directory already exists: %s\n' % directory_path)
	else:
		print ('Making new directory: %s\n' % directory_path)
		os.mkdir(directory_path)

# Get a list of the accessions
def get_accessions(directory):

    accessions = []
    for file in os.listdir(directory):
        if not file.endswith('.DS_Store'):
            accessions.append(re.findall(r'(.+)(?=\.)', file)[0])

    return accessions

# Get a list of the cds sequences for a genome
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


# Parrallelise
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

# Get the off frame freqeuncies of each codon
def get_off_frame_frequencies(seqs):

    off_frame_frequencies = {}

    codon_count = 0

    for frame in frames:
        off_frame_frequencies[frame] = {}
        for codon in codon_list:
            off_frame_frequencies[frame][codon] = 0

    for cds in seqs:
        codon_count += (len(cds)-3)/3

        f1 = re.findall(r'.{3}', cds[1:])
        f2 = re.findall(r'.{3}', cds[2:])
        for codon in codon_list:
            off_frame_frequencies[1][codon] += f1.count(codon)
            off_frame_frequencies[2][codon] += f2.count(codon)

    return(off_frame_frequencies, codon_count)

# Calculate the GC content of the sequences
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


# Get the genome stats
def genome_stats(cds_seqs):

    gc, gc3 = calc_gc(cds_seqs)
    off_frame_counts, codon_count = get_off_frame_frequencies(cds_seqs)
    return(gc, gc3, off_frame_counts, codon_count)


# Get the codons in the sequences
def get_cds_codons(cds_seqs):

    cds_codons = []

    for cds in cds_seqs:
        codons = re.findall(r'.{3}', cds)
        cds_codons.append(codons)
    return(cds_codons)

# Generate shuffled sequences
def shuffle_genome_codons(cds_codons):

    shuffled_seqs = []

    for cds in cds_codons:

        # Shuffle the codons, excluding the start and stop
        shuffled = shuffle_codons(cds[1:-1])

        # Create a string of the new sequence
        new_seq = cds[0] + conconate(shuffled, '') + cds[-1]

        # Add the new sequence to the array
        shuffled_seqs.append(new_seq)

    return(shuffled_seqs)



def scramble(a, axis=-1):
    """
    Return an array with the values of `a` independently shuffled along the
    given axis
    """
    b = np.random.random(a.shape)
    idx = np.argsort(b, axis=axis)
    shuffled = a[np.arange(a.shape[0])[:, None], idx]
    return (shuffled)


def join(codons):

    joined =''.join(codons)
    return(joined)

# Run the randomisation for the required number of repeats
def run_repeats(repeat_list, cds_codons, accession):

    outputs = []

    # Set of codons without start and stop codons
    without = [codon_set[1:-1] for codon_set in cds_codons]

    # Get max number of codons in the cds
    top = max(len(codon_set) for codon_set in without)

    # Generate array of start and stop codons and transform
    starts = np.array([codon_set[0] for codon_set in cds_codons])
    stops = np.array([codon_set[-1] for codon_set in cds_codons])

    starts = np.reshape(starts, (len(starts),1))
    stops = np.reshape(stops, (len(stops),1))

    # Set up codon set with blanks to get correct dimensions (each cds)
    # results in the same number of 'codons', but filled with blanks
    spaces = [codon_set + (['']*(top-len(codon_set))) for codon_set in without]
    spaces = np.array(spaces)


    for repeat in repeat_list:

        np.random.seed()

        print('Repeat %s' % (repeat+1))

        # Shuffle the codons
        t = scramble(spaces)

        # Merge the start and stop codons back with the shuffled sequences
        merge = np.hstack((starts, t, stops))

        # Concatonate the sequences
        joined = [join(codons) for codons in merge]

        # Get the stats for the new sequences
        gc, gc3, off_frame_counts, codon_count = genome_stats(joined)

        # Return the outputs
        output = [accession, repeat+1, gc, gc3, off_frame_counts, codon_count]
        outputs.append(output)

    return(outputs)


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
def write_repeats(output_file, results):

    # For each of the repeats, write the line
    for result in results:

        outputs = result.get()

        for output in outputs:

            acc = output[0]
            repeat = output[1]
            gc = output[2]
            gc3 = output[3]
            off_frame_counts = output[4]
            codon_count = output[5]

            line = '%s_r%s,%s,%s,%s' % (acc, repeat, gc, gc3, codon_count)

            for frame in sorted(off_frame_counts):
                for codon in sorted(off_frame_counts[frame]):
                    line += ',%s' % ((100/codon_count)*off_frame_counts[frame][codon])

            line += '\n'
            output_file.write(line)


# Write the results to file
def write_to_file(accession, real_stats, results):

    output_file = open('outputs/simulation_codon_shuffle/' + accession + '_codon_shuffle.csv', 'w')

    write_header(output_file)
    write_real(output_file, accession, real_stats)
    write_repeats(output_file, results)

    output_file.close()



def run_genomes(accessions, acc_counts):


    for accession in accessions:

        print('%s: %s' % (acc_counts[accession], accession))

        # Get the sequences for the genome
        cds_seqs = get_seqs(accession)

        # cds_seqs = cds_seqs[:1]

        cds_codons = get_cds_codons(cds_seqs)


        # Get gc, gc3, off frame counts for the real sequences
        gc, gc3, real_off_counts, codon_count = genome_stats(cds_seqs)
        real_stats = [gc, gc3, real_off_counts, codon_count]

        repeat_list = list(range(required_repeats))
        results = run_in_parralell(repeat_list, [cds_codons, accession], run_repeats)
        # run_repeats(repeat_list, cds_codons, accession)

        # Write the results to file
        write_to_file(accession, real_stats, results)



def main():


    t0_main = time.time()

    create_directory('outputs/simulation_codon_shuffle/')

    # Return the accessions of the genomes in the directory
    accessions = get_accessions('outputs/genome_extractions/')

    genomes = []
    acc_counts = {}


    accession_count = 0

    for accession in sorted(accessions):
        accession_count += 1
        if accession_count:
        # if accession_count <= 1:
        # if accession == 'CP003881':
            genomes.append(accession)
            acc_counts[accession] = accession_count

    print('\nShuffling codons within the CDS for %s genomes\n' % len(genomes))

    # results = run_in_parralell(genomes, [acc_counts], run_genomes)
    run_genomes(genomes, acc_counts)


    t1_main = time.time()
    print ('\n%s\nTotal time: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))


if __name__ == '__main__':
    main()
