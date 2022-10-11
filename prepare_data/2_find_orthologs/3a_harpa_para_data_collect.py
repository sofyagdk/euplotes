####################################
# sofya 08.04.2022
# sofya.gaydukova@gmail.com

# read aligned files from ./harpa_para_1turned/
# this folder is alreay manually cleaned from 
# harpa orthogroups containing more than two sequencies 
# > 2  are stored in ./harpa_many
#
# calculate identity
# create table with identity and lengths
#
# results are printed into STDOUT
# store in > para_stats_v2.txt
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


def CalcATComposition(seq):
	seq = seq.replace('-', '')
	at = 0
	for letter in seq:
		if letter in ['A', 'T']:
			at += 1
	return at, len(seq)



def CalcMultipleIdentity(alignment_dict):
	alignment = list()
	for idr in alignment_dict:
		alignment.append(list(alignment_dict[idr]))
	alignment = np.array(alignment).T
	# alignment = alignment.T

	conservative = 0 
	non_gaped = 0
	
	for col in alignment:
		if '-' in col:
			continue
		non_gaped += 1
		if len(set(col)) == 1:
			conservative += 1
	return conservative, non_gaped





####################################
# DATA: tables and dictionaries
####################################

alignment_dir = './harpa_para_1turned/'

alignment_files = os.listdir(alignment_dir )


print('filename\tlen\tcols_conservative\tcols_nongap\tsp_idrs\tsp_lens')
	


for filename in alignment_files:
	path = alignment_dir+filename
	bp_dict = reader.readfasta_todictionary(path)
	conservative, non_gaped = CalcMultipleIdentity(bp_dict)

	idrs = bp_dict.keys()
	lens = []

	for idr in idrs:
		lens.append(str(len(bp_dict[idr].replace('-', ''))))


	print('\t'.join([filename, str(len(bp_dict)), str(conservative), str(non_gaped), ','.join(idrs), ','.join(lens)]))