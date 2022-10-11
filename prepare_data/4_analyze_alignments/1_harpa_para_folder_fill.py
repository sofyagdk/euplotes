####################################
# sofya 08.04.2022
# sofya.gaydukova@gmail.com

# make HARPA paralog files
# from proteinortho result file 
# uses harpa_para_table.csv
# created files are added to ./harpa_para/
# files are not aligned
# next file -- para_turn.py
# that turns and alignes files into folder ./harpa_para_1turned/
####################################


import os
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader
import pandas as pd 
import numpy as np  




####################################
# DATA: tables and dictionaries
####################################


transcriptome_path = '/mnt/gamma/user/sofya/data/transcriptomes/filtered/fE.harpa_n.fasta'


transcriptome_dict = reader.readfasta_todictionary(transcriptome_path)
output_folder = './harpa_para/'



table = pd.read_csv('harpa_para_table.csv', header = 0, sep = '\t')
####################################
# MAIN part: read and write
####################################



for para_pair_line in table['fE.harpa_n.fasta']:
	para_pair = para_pair_line.split(',')
	with open(output_folder + para_pair[0] + '.fasta', 'w') as fout:
		for para in para_pair:
			fout.write(''.join(['>', para, '\n'] ))
			fout.write(transcriptome_dict[para] + '\n')
