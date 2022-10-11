####################################
# sofya 08.04.2022
# sofya.gaydukova@gmail.com

# 
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

alignment_dir = '../nuc_alignments/'

alignment_files = os.listdir(alignment_dir )


# print('filename\tlen\tcols_conservative\tcols_nongap\tsp_idrs\tsp_ats\tsp_nongap\tsp_atcomp')
	

# for filename in alignment_files:
# 	path = alignment_dir+filename
# 	bp_dict = reader.readfasta_todictionary(path)
# 	conservative, non_gaped = CalcMultipleIdentity(bp_dict)

# 	idrs = ','.join(bp_dict.keys())
# 	ats = list()
# 	nongaps = list()
# 	atcomp = list()
# 	for idr in bp_dict.keys():
# 		at, nongap = CalcATComposition(bp_dict[idr])
# 		ats.append(str(at))
# 		nongaps.append(str(nongap))
# 		atcomp.append(str(at/float(nongap)))



# 	print('\t'.join([filename, str(len(bp_dict)), str(conservative), str(non_gaped), idrs, ','.join(ats), ','.join(nongaps), ','.join(atcomp)]))




print('filename\tlen\tcols_conservative\tcols_nongap\tsp_idrs\tsp_ats\tsp_nongap\tsp_atcomp')
	

for filename in alignment_files:
	path = alignment_dir+filename
	bp_dict = reader.readfasta_todictionary(path)
	idrs = bp_dict.keys()

	for i_key in range(len(idrs)):
		for j_key in range(i_key+1, len(idrs)):

			small_dict = dict()
			small_dict[idrs[i_key]] = bp_dict[idrs[i_key]]
			small_dict[idrs[j_key]] = bp_dict[idrs[j_key]]

			conservative, non_gaped = CalcMultipleIdentity(small_dict)

			small_idrs = ','.join(small_dict.keys())
			ats = list()
			nongaps = list()
			atcomp = list()
			
			for idr in small_dict.keys():
				at, nongap = CalcATComposition(small_dict[idr])
				ats.append(str(at))
				nongaps.append(str(nongap))
				atcomp.append(str(at/float(nongap)))



			print('\t'.join([filename, str(len(small_dict)), str(conservative), str(non_gaped), small_idrs,  ','.join(ats), ','.join(nongaps), ','.join(atcomp)]))
