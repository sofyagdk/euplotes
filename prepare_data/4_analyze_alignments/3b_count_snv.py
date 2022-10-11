####################################
# sofya 11.04.2022
# sofya.gaydukova@gmail.com

# calculate SNVs
# prints into STDOUT
# > ./snv_data/sp1_sp2_snvs.txt
####################################



import os
import math
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader
import numpy as np

####################################
# FUNCTIONS
####################################


def FindSNVs(alignment_dict, ):

	alignment = list()
	for idr in alignment_dict:
		alignment.append(list(alignment_dict[idr]))
	alignment = np.array(alignment).T

	snvs = list()

	
	for index, col in enumerate(alignment):
		if '-' in col:
			continue
		if col[0] != col[1]:
			snvs.append(list(col) + [str(index)])

	return snvs

####################################
# MAIN -- for HARPA
####################################





alignment_dir = './harpa_para_1t/'
filenames = os.listdir(alignment_dir)


# for filename in filenames:
# 	path = alignment_dir + filename
# 	bp_dict = reader.readfasta_todictionary(path)
# 	snvs = FindSNVs(bp_dict)
# 	for snv in snvs:
# 		print(' '.join([filename] + list(snv)))


####################################
# MAIN -- for other two species
####################################




alignment_dir = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/nuc_gold/'

filenames = os.listdir(alignment_dir)



for filename in filenames:
	path = alignment_dir + filename
	bp_dict = reader.readfasta_todictionary(path)

	interesting = dict()
	for key in bp_dict:
		if 'eury' in key:
			interesting[key] = bp_dict[key]
		if 'rari' in key:
			interesting[key] = bp_dict[key]
	if len(interesting) < 2:
		continue


	snvs = FindSNVs(interesting)
	for snv in snvs:
		print(' '.join([filename] + list(snv)))
