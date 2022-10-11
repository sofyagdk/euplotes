
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import evo_ans_temp_v4 as gardener
from ete3 import Tree
import reader
import os 
import subprocess
import bootstrap_tables_v4 as boot 
import re
import numpy as np


codingsfilename = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/mark.txt'

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


#############
## INSERTED
#############

codon_list = [
"TTT", "TCT", "TAT", "TGT", 
"TTC", "TCC", "TAC", "TGC", 
"TTA", "TCA", "TAA", "TGA", 
"TTG", "TCG", "TAG", "TGG", 
"CTT", "CCT", "CAT", "CGT", 
"CTC", "CCC", "CAC", "CGC", 
"CTA", "CCA", "CAA", "CGA", 
"CTG", "CCG", "CAG", "CGG", 
"ATT", "ACT", "AAT", "AGT", 
"ATC", "ACC", "AAC", "AGC", 
"ATA", "ACA", "AAA", "AGA", 
"ATG", "ACG", "AAG", "AGG", 
"GTT", "GCT", "GAT", "GGT", 
"GTC", "GCC", "GAC", "GGC", 
"GTA", "GCA", "GAA", "GGA", 
"GTG", "GCG", "GAG", "GGG" 
]


# table:
# codon next_letter such_events such_codon_in_sp species

# dict- 
#  sp:codon_dict
# codon_dict:
#  codon: next_letter_list


codon_dict = {'petz':dict(),
'raik':dict(),
'octo':dict(),
'harp':dict(),
'eury':dict(),
'rari':dict(),
'foca':dict(),
'minu':dict(),
'cras':dict()}

species_list = ['petz','raik','octo','harp','eury','rari','foca','minu','cras']

for codon in codon_list:
	for sp in species_list:
		codon_dict[sp][codon] = list()



directory ='/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/nuc_gold/'
files = os.listdir(directory)



for afile in files:
	filename = directory + afile
	before_dict = reader.readfasta_todictionary(filename)

	for idr in before_dict:
		before_dict[idr]= boot.clean(before_dict[idr])

	for idr in before_dict:
		animal = idr[0:4]


		amount = 0
		coding = fish_coding_true(idr)
		seq = boot.clean(before_dict[idr])
		if coding == list():
			continue

		

		for i in range(0,len(coding),2):
			translated_piece = seq[coding[i]:coding[i + 1]]

				
			for i in range(0, len(translated_piece) - 6, 3):
				curr_codon = translated_piece[i:i+3]
				codon_dict[animal][curr_codon].append(translated_piece[i+3])
				


# codon next_letter such_events such_codon_in_sp species



with open('AAA_codon_stats.txt', 'a') as fout:
	fout.write('animal codon next_letter such_events num_such_codon_in_sp')
	for sp in species_list:
		for codon in codon_list:
			for letter in ['A', 'T', 'G', 'C']:

				such_events = str(np.sum(np.array(codon_dict[sp][codon]) == letter))
				num_such_codon_in_sp = str(len(codon_dict[sp][codon]))

				s = ' '.join([sp, codon, letter, such_events, num_such_codon_in_sp ]) 
				s += '\n'
				fout.write(s)



