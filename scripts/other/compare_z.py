#!/usr/bin/python

# Script number:
# File:
# Prerequisite script(s):		simulation scripts
# Prerequisite file(s):
# Description:					Comparse z scores
# Output files:                 compare_z.csv

import numpy as np
import time
import os
import rpy2.robjects as robjects

def create_directory(directory_path):

	if os.path.exists(directory_path):
		print ('Directory already exists: %s\n' % directory_path)
	else:
		print ('Making new directory: %s\n' % directory_path)
		os.mkdir(directory_path)

def array_to_arraystring(array):
    string = 'c('
    for i in array:
        string += '{},'.format(i)
    string = string[:-1] + ')'
    return(string)





def get_z_scores(path):

	stops = ['TAA', 'TAG', 'TGA']

	p1_zs = {}
	p1_ps = {}

	for stop in stops:
		p1_zs[stop] = {}
		p1_ps[stop] = []

	accs = []

	with open(path, 'rU') as my_file:


		lines = my_file.readlines()

		for line in lines[1:]:



			splits = line.split(',')
			acc = splits[0]
			accs.append(acc)


			TAA1 = float(splits[3])
			TAG1 = float(splits[5])
			TGA1 = float(splits[7])

			p1_zs['TAA'][acc] = TAA1
			p1_zs['TAG'][acc] = TAG1
			p1_zs['TGA'][acc] = TGA1

			p1_ps['TAA'].append(float(splits[4]))
			p1_ps['TAG'].append(float(splits[6]))
			p1_ps['TGA'].append(float(splits[8]))


	p1_padj = {}
	for stop in stops:
		p1_padj[stop] = {}

	for stop in p1_ps:
		padj = robjects.r('p.adjust({}, method="fdr")'.format(array_to_arraystring(p1_ps[stop])))

		for i in range(len(padj)):
			p1_padj[stop][accs[i]] = padj[i]


	return(p1_zs, p1_padj)

	# with open('outputs/simulation_codon_shuffle_analysis/stops.csv', 'rU') as my_file:
    #
	# 	lines = my_file.readlines()
    #
	# 	for line in lines[1:]:
    #
	# 		splits = line.split(',')
	# 		acc = splits[0]
    #
	# 		TAA1 = float(splits[3])
	# 		TAG1 = float(splits[5])
	# 		TGA1 = float(splits[7])
    #
	# 		p1_zs['TAA'][acc].append(TAA1)
	# 		p1_zs['TAG'][acc].append(TAG1)
	# 		p1_zs['TGA'][acc].append(TGA1)
    #
	# return(p1_zs)



def compare(p1_zs_ss, p1_padj_ss, p1_zs_cs, p1_padj_cs):

	output_file = open('outputs/other/compare_z.csv', 'w')
	output_file.write('codon,sig_pos_codon_shuffle,sig_pos_syn_site,sig_pos_both,%_both/max(cs/ss)\n')

	for codon in p1_zs_ss:

		positive_ss = []
		positive_cs = []
		same = []

		for acc in p1_zs_ss[codon]:
			if p1_zs_ss[codon][acc] > 0 and p1_padj_ss[codon][acc] < 0.05:
				positive_ss.append(acc)
			if p1_zs_cs[codon][acc] > 0 and p1_padj_cs[codon][acc] < 0.05:
				positive_cs.append(acc)
			if p1_zs_ss[codon][acc] > 0 and p1_padj_ss[codon][acc] < 0.05 and p1_zs_cs[codon][acc] > 0 and p1_padj_cs[codon][acc] < 0.05:
				same.append(acc)

		output_file.write('{},{},{},{},{}\n'.format(codon, len(positive_cs), len(positive_ss), len(same), np.divide(len(same), max(len(positive_cs), len(positive_ss)))))


	output_file.close()

def run():

	p1_zs_ss, p1_padj_ss = get_z_scores('outputs/simulation_synonymous_site_analysis/stops.csv')
	p1_zs_cs, p1_padj_cs = get_z_scores('outputs/simulation_codon_shuffle_analysis/stops.csv')

	compare(p1_zs_ss, p1_padj_ss, p1_zs_cs, p1_padj_cs)



def main():

	t0_main = time.time()

	create_directory('outputs/other/')
	run()

	t1_main = time.time()
	print ('\n%s\nTotal time: %s\n%s\n' % ('='*30, t1_main-t0_main, '='*30))


if __name__ == '__main__':
	main()
