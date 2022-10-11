# last change 20-08-04 to count the contexts in animals into file for_statistics.txt

import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import evo_ans_temp_v4 as gardener
from ete3 import Tree
import reader
import os 
import subprocess
# import marking_module_v3 as mark
import bootstrap_tables_v4 as boot 
import re


codingsfilename = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/mark.txt'
# codingsfilename = '/mnt/gamma/user/sofya/scripts/final/0pipeline/18_graphmark/mark_raw.txt'

def fish_coding_alignment(idr):
	global codingsfilename
	with open(codingsfilename, 'r') as fin:
		for line in fin:
			linelist = line.split()
			# cras_c3120_g1_i1	[35, 842, 843, 885]	[40, 848, 849, 920]	[40, 848]
			if linelist[1] == idr:
				cod = linelist[4]
				cod = list(map(int, cod.split(',')))
				return cod

def fish_coding_true(idr):
	global codingsfilename
	with open(codingsfilename, 'r') as fin:
		for line in fin:
			linelist = line.split()
			# cras_c3120_g1_i1	[35, 842, 843, 885]	[40, 848, 849, 920]	[40, 848]
			if linelist[1] == idr:
				cod = linelist[2]
				cod = list(map(int, cod.split(',')))
				return cod


# #############
# ## INSERTED
# #############

# cont_dict = {'petz':[0, 0, 0, 0, 0, 0, 0],
# 'raik':[0, 0, 0, 0, 0, 0, 0],
# 'octo':[0, 0, 0, 0, 0, 0, 0],
# 'harp':[0, 0, 0, 0, 0, 0, 0],
# 'eury':[0, 0, 0, 0, 0, 0, 0],
# 'rari':[0, 0, 0, 0, 0, 0, 0],
# 'foca':[0, 0, 0, 0, 0, 0, 0],
# 'minu':[0, 0, 0, 0, 0, 0, 0],
# 'cras':[0, 0, 0, 0, 0, 0, 0]}



# directory ='/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/nuc_gold/'
# files = os.listdir(directory)


# # contexts_stop = ['AA', 'AG']

# for afile in files:
# 	filename = directory + afile
# 	before_dict = reader.readfasta_todictionary(filename)

# 	for idr in before_dict:
# 		before_dict[idr]= boot.clean(before_dict[idr])

# 	for idr in before_dict:
# 		animal = idr[0:4]


# 		amount = 0
# 		coding = fish_coding_true(idr)
# 		seq = boot.clean(before_dict[idr])
# 		if coding == list():
# 			continue

		
# 		cont_dict[animal][0] += len(coding)//2 - 1


# 		search_region = str()
# 		pruned_search_region = str()
# 		for i in range(0,len(coding),2):
# 			search_region += seq[coding[i]:coding[i + 1]]
# 			pruned_search_region +=  seq[coding[i] + 15 : coding[i + 1] - 15]
		
# 		cont_dict[animal][6] += len(pruned_search_region)
# 		cont_dict[animal][5] += len(search_region)

				
# 		for i in range(0, len(search_region) - 6, 3):
# 			if search_region[i:i+2] == 'AA':
# 				cont_dict[animal][1] += 1
# 			elif search_region[i:i+2] == 'AG':
# 				cont_dict[animal][3] += 1


# 		for i in range(0, len(pruned_search_region) - 6, 3):
# 			if search_region[i:i+2] == 'AA':
# 				cont_dict[animal][2] += 1
# 			elif search_region[i:i+2] == 'AG':
# 				cont_dict[animal][4] += 1


# print(cont_dict)
# with open('animal_statistics_stops.txt', 'a') as fout:
# 	fout.write('animal\tshifts\tAA\tAA_15\tAG\tAG_15\tsum_length\tsum_length_15\n')
# 	for idr in cont_dict:
# 		s = '\t'.join(list(map(str, [idr, cont_dict[idr][0], 
# 			cont_dict[idr][1], cont_dict[idr][2], cont_dict[idr][3], 
# 			cont_dict[idr][4], cont_dict[idr][5], cont_dict[idr][6]])))
# 		s += '\n'
# 		fout.write(s)






# exit()
# #############
# ## END
# #############



cont_dict = {'petz':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
'raik':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
'octo':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
'harp':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
'eury':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
'rari':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
'foca':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
'minu':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
'cras':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]}

# directory ='/mnt/gamma/user/sofya/scripts/final/homo_bp_al_3_v8/'
# directory ='/mnt/gamma/user/sofya/scripts/final/homo_bp_al_3_v2_0/'
# directory ='/mnt/gamma/user/sofya/scripts/final/homo_bp_rev_alignedwater/'
directory ='/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/nuc_gold/'
# directory ='/mnt/gamma/user/sofya/scripts/final/0pipeline/nuc_alignments/'
files = os.listdir(directory)

# contexts_1 = ['AA', 'GA']
# contexts_2 = ['AA', 'GA', ]
contexts_AA = ['AAAAA','AAAAG']
contexts_AG = ['AAGAA','AAGAG']


for afile in files:
	filename = directory + afile
	before_dict = reader.readfasta_todictionary(filename)

	for idr in before_dict:
		before_dict[idr]= boot.clean(before_dict[idr])

	# sp_list = list()
	# paralouge_list = list()

	# for idr in before_dict:
	# 	animal = idr[0:4]
	# 	if animal in sp_list:
	# 		paralouge_list.append(animal)
	# 	sp_list.append(animal)

	for idr in before_dict:
		animal = idr[0:4]
		# if animal in paralouge_list:
		# 	continue
		amount = 0
		coding = fish_coding_true(idr)
		seq = boot.clean(before_dict[idr])
		if coding == list():
			continue

		
		cont_dict[animal][0] += len(coding)//2 - 1
		cont_dict[animal][6] += 1


###############################################################

		for i in range(2, len(coding), 2):
			c = seq[coding[i-1] - 3:coding[i-1]]
			if c == 'AAA':
				cont_dict[animal][7] += 1
			if c == 'AAG':
				cont_dict[animal][8] += 1


###############################################################


		search_region = str()
		for i in range(0,len(coding),2):
			search_region += seq[coding[i]:coding[i + 1]]
		
		cont_dict[animal][5] += len(search_region)

				
		for i in range(0, len(search_region) - 6, 3):
			if search_region[i:i+5] in contexts_AA:
				cont_dict[animal][1] += 1
			elif search_region[i:i+5] in contexts_AG:
				cont_dict[animal][3] += 1

		for i in range(0, len(search_region) - 6):
			if search_region[i:i+5] in contexts_AA:
				cont_dict[animal][2] += 1
			elif search_region[i:i+5] in contexts_AG:
				cont_dict[animal][4] += 1

		for i in range(0, len(search_region) - 6, 3):
			if search_region[i:i+2] == 'AA':
				cont_dict[animal][9] += 1
			elif search_region[i:i+2] == 'AG':
				cont_dict[animal][10] += 1


print(cont_dict)
with open('./animal_statistics_gold_long.txt', 'a') as fout:
	fout.write('animal  shifts  AA_or  AA_no  AG_or  AG_no  sum_length sum_transcripts per_codon*100000 per_transript*100 AAAshifts AAGshifts AA AG\n')
	for idr in cont_dict:
		s = '\t'.join(list(map(str, [idr, cont_dict[idr][0], cont_dict[idr][1], 
			cont_dict[idr][2], cont_dict[idr][3], cont_dict[idr][4], cont_dict[idr][5], 
			cont_dict[idr][6], cont_dict[idr][0]/float(cont_dict[idr][5]/3)*100000, cont_dict[idr][0]/float(cont_dict[idr][6])*100, 
			cont_dict[idr][7], cont_dict[idr][8], 
			cont_dict[idr][9], cont_dict[idr][10]])))
		s += '\n'
		fout.write(s)





# t = Tree('(petz,(raik,((octo,harp),(eury,(rari,(foca,(minu,cras)))))));')