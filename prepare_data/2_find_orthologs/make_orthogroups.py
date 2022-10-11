####################################
# sofya 08.04.2022
# sofya.gaydukova@gmail.com

# make orthogroup files
# from proteinortho result file 
# filters paralogues
# also include harpa_to_take.txt
# that is created with /home/sofya/main/Lab/Euplotes/16_protein_ortho/protortho_analyse.ipynb
# (just the longest seq from two paralogous harpas)
####################################


import os
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader
import numpy as np
import pandas as pd



####################################
# DATA: tables and dictionaries
####################################


transcriptome_directory = '/mnt/gamma/user/sofya/data/transcriptomes/filtered/'


transcriptome_fastas = ['fE.euryhalinus_n.fasta', 
'fE.focardii_n.fasta', 
'fE.harpa_n.fasta', 
'fE.minuta_n.fasta', 
'fE.octocarinatus_n.fasta', 
'fE.petzi_n.fasta', 
'fE.raikovi_n.fasta', 
'fE.rariseta_n.fasta', 
'fT.thermophila_n.fasta', 
'fE.crassus_combined.fasta']


transcriptome_dict = dict()

for tr in transcriptome_fastas:
	if 'thermophi' in tr:
		continue
	path = transcriptome_directory + tr
	species = tr[3:7]
	transcriptome_dict[species] = reader.readfasta_todictionary(path)




path_to_proteinortho_file = '/mnt/gamma/user/sofya/proteinortho_results/filtered_e25_id70.proteinortho'
output_folder = '../orthogroups/'

####################################
# MAIN part: read and write
####################################


# Species	Genes	Alg.-Conn.	fE.crassus_combined.fasta	fE.euryhalinus_n.fasta	fE.focardii_n.fasta	fE.harpa_n.fasta	fE.minuta_n.fasta	fE.octocarinatus_n.fasta	fE.petzi_n.fasta	fE.raikovi_n.fasta	fE.rariseta_n.fasta	fT.thermophila_n.fasta

# 7	8	0.125	cras_c3754_g1_i1++	eury_eury_c19129_g1_i1	foca_c32256_g1_i1	harp_c1029_g1_i1,harp_c1457_g1_i1	minu_c34451_g1_i1	octo_c70426_g1_i1	*	*	rari_c876_g1_i1	*



harpa_to_take_table = pd.read_csv('harpa_to_take.txt', header = 0, sep=' ')
harpa_to_take = list(harpa_to_take_table['to_take'])








with open(path_to_proteinortho_file, 'r') as fin:
	for line in fin:
		if line[0] == '#':
			continue

		orthogroup = list()

		line_split = line.strip().split('\t')
		for idr in line_split[3:]:
			if ',' in idr: ######### means paralog 
				
				para_idrs = idr.split(',')
				for para_idr in para_idrs:
					if para_idr in harpa_to_take:
						orthogroup.append(para_idr)
						print('yo!')
				
				continue
			
			if idr == '*':
				continue
			
			orthogroup.append(idr)

		

		if len(orthogroup) < 2:
			print('too small', orthogroup)
			continue
		
		
		with open(output_folder + orthogroup[0] + '.fasta', 'w') as fout:
			for idr in orthogroup:
				sp = idr[:4]
				fout.write(''.join(['>', idr, '\n']))
				fout.write(transcriptome_dict[sp][idr])
				fout.write('\n')









