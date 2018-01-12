#!/usr/bin/python

# Script number:				7.9
# File:						    9 of 9
# Prerequisite script(s):       synonymous_codon_simulation, synonymous_codon_simulation_analysis
# Prerequisite file(s):
# Description:					Output excess stats of model
# Output files:                 stats.csv



import numpy as np
import re
import time
import multiprocessing as mp
import os
import sys
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()



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

frames = [1, 2, 'both']

def read_combined_stops():

    zs = {}
    ps = {}
    gc = []
    gc3 =[]

    for frame in frames:
        zs[frame] = []
        ps[frame] = []

    with open('outputs/simulation_synonymous_codon_analysis/combined_stops.csv', 'rU') as file:
        lines = file.readlines()

        for line in lines[1:]:
            line = line.strip('\n').split(',')
            zs[1].append(float(line[3]))
            zs[2].append(float(line[6]))
            zs['both'].append(float(line[9]))
            ps[1].append(float(line[5]))
            ps[2].append(float(line[8]))
            ps['both'].append(float(line[11]))
            gc.append(float(line[1]))
            gc3.append(float(line[2]))

    gc = array_to_arraystring(gc)
    gcs_string = array_to_arraystring(gc)

    outputs = {}

    for frame in frames:

        pos = calc_pos_z(zs[frame], ps[frame])
        cor = calc_cor(gc, zs[frame])
        outputs[frame] = [pos, len(zs[frame]), cor]

    return(outputs)

def calc_pos_z(zs, ps):

    padj = robjects.r('p.adjust({}, method="fdr")'.format(array_to_arraystring(ps)))

    pos = 0
    for i in range(len(zs)):
        if(zs[i] > 0 and padj[i] < 0.05):
            pos += 1
    return(pos)

def calc_cor(gc, zs):

    cor = robjects.r('cor.test({}, {}, method="spearman")'.format(gc, array_to_arraystring(zs)))
    return(extractFloatVetor(cor[3]), extractFloatVetor(cor[2]))

def read_codon_outputs():

    query_codons = ['TAA', 'TAC', 'TAG', 'TAT', 'TGA', 'TGC', 'TGG', 'TGT']
    # query_codons = ['TAG']

    zs = {}
    ps = {}
    gc = []
    gc3 =[]

    for frame in frames:
        zs[frame] = {}
        ps[frame] = {}

        for codon in query_codons:
            zs[frame][codon] = []
            ps[frame][codon] = []

    with open('outputs/simulation_synonymous_codon_analysis/all_codons.csv', 'rU') as file:
        lines = file.readlines()

        for line in lines[1:]:
            line = line.strip('\n').split(',')
            gc.append(float(line[1]))
            gc3.append(float(line[2]))

            for codon in query_codons:
                index_pos = codon_list.index(codon)
                index = 2*index_pos

                zs[1][codon].append(float(line[index+3]))
                ps[1][codon].append(float(line[index+4]))
                zs[2][codon].append(float(line[index+131]))
                ps[2][codon].append(float(line[index+132]))
                zs['both'][codon].append(float(line[index+259]))
                ps['both'][codon].append(float(line[index+260]))

    gc = array_to_arraystring(gc)
    gc3 = array_to_arraystring(gc)


    outputs = {}
    for frame in frames:
        outputs[frame] = {}

    for frame in frames:
        for codon in sorted(zs[frame]):
            pos = calc_pos_z(zs[frame][codon], ps[frame][codon])
            cor = calc_cor(gc, zs[frame][codon])
            outputs[frame][codon] = [pos, len(zs[frame][codon]), cor]

    return(outputs)



def extractFloatVetor(item):
    return(item[0])

def array_to_arraystring(array):
    string = 'c('
    for i in array:
        string += '{},'.format(i)
    string = string[:-1] + ')'
    return(string)


def write_to_file(combined, codons):

    output_file = open('outputs/simulation_synonymous_codon_analysis/stats.csv', 'w')
    output_file.write('codon,frame,num_excess,prop_with_excess,cor_rho,cor_p\n')

    for frame in combined:
        output = combined[frame]
        output_file.write('all_stops,{},{},{},{},{}\n'.format(frame, output[0], np.divide(output[0], output[1]), output[2][0], output[2][1]))

    for frame in codons:
        for codon in codons[frame]:
            output = codons[frame][codon]
            output_file.write('{},{},{},{},{},{}\n'.format(codon, frame, output[0], np.divide(output[0], output[1]), output[2][0], output[2][1]))

    output_file.close()


def main():

    t0_main = time.time()

    combined_outputs = read_combined_stops()
    codon_outputs = read_codon_outputs()
    write_to_file(combined_outputs, codon_outputs)




    t1_main = time.time()
    print ('\n%s\nTotal time: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))


if __name__ == '__main__':
    main()

# combined_stops_excess <- function(frame) {
#   file <- read.csv('outputs/simulation_codon_shuffle_analysis/combined_stops.csv', head=T)
#
#   col <- paste('osc_', frame, '_z', sep='')
#   pcol <- paste('osc_', frame, '_pval', sep='')
#
#   file$padj <- p.adjust(file[[pcol]], method="fdr")
#
#   print(cor.test(file$gc, file[[col]], method="spearman"))
#   print(nrow(file[file[[col]] > 0 & file$padj < 0.05,]))
#   print(nrow(file[file[[col]] > 0 & file$padj < 0.05,])/nrow(file)*100)
#
# }
