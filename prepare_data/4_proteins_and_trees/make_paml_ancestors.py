# this needed for the algorithmic part
import random
import os
import glob
import subprocess
import numpy as np

# this needed for the main part
import os
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader
import random 
from optparse import OptionParser
import ds




#####################################################
### Dictionaries and Data
#####################################################

codings_true = dict()

with open('/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/mark.txt', 'r') as fin:
	for line in fin:
		cod = line.split()[2]
		cod = list(map(int, cod.split(',')))
		codings_true[line.split()[1]] = cod



#####################################################
### Functions
#####################################################
def create_filename(flag = None, extention = '.txt', num = 15):
	filename = str(flag)
	for i in range(num):
		filename += random.choice('0123456789') 
	filename += extention
	return filename



def get_paml_tree(codon_filename, tree_filename):

	dirname =  create_filename(flag = "paml_temp", extention = '', num = 15)
	os.mkdir('./' + dirname)
	os.system('cp /mnt/gamma/user/sofya/scripts/final/0pipeline/13_LysEvo/codon_alignments_nuc/'+ codon_filename +' ./' + dirname + '/codon_alignment.nuc')
	os.system('cp /mnt/gamma/user/sofya/scripts/final/0pipeline/13_LysEvo/nwk_trees/'+ tree_filename +' ./' + dirname + '/tree_file.nwk')
	os.system('cp ./codeml.clt ./' + dirname + '/codeml.clt')

	# os.system("sed -i 's/outfilename/.\/paml_ancestors\/" + codon_filename + "/g' " + dirname + '/codeml.clt')
	os.system("sed -i 's/codon_alignment.nuc/" + dirname + "\/codon_alignment.nuc/g' "+ dirname + '/codeml.clt')
	os.system("sed -i 's/tree_file.nwk/" + dirname + "\/tree_file.nwk/g' "+ dirname + '/codeml.clt')
	# seqfile = ./two_seqs.nuc * sequence data file name 
	# outfile = codeml_2_seqs * main result file)

	# pause = int(input('Check directory ----> '))



	os.system('codeml /mnt/gamma/user/sofya/scripts/final/0pipeline/13_LysEvo/' + dirname + '/codeml.clt') #  > ./' + dirname + '/codeml.out')

	
	# pause = int(input('Paml has finished ----> '))
	animals = ['harp', 'petz', 'octo', 'minu', 'cras', 'foca', 'raik', 'rari', 'eury']
	outdict = dict()
	
	with open('rst', 'r') as fin:
		read_tree_info = False
		read_seqs = False
		for line in fin:
			if line == "tree with node labels for Rod Page's TreeView\n":
				read_tree_info = True 
				continue
			elif read_tree_info:
				tree_info = line.strip()
				read_tree_info = False
				continue
			elif line == "List of extant and reconstructed sequences\n":
				read_seqs = True
				continue
			elif read_seqs:
				if line[0] == '\n' or line[0] == ' ':
					continue
				if "Overall accuracy of the" in line:
					read_seqs = False
					break

				if line[:4] in animals:
					outdict[line[:4]] = line[4:-1].replace(" ", "")
				elif line[:4] == 'node':
					outdict[line[:9].replace(" ", "")] = line[9:-1].replace(" ", "")



	os.system('rm -r '+ dirname)
	return tree_info, outdict





#####################################################
### Main part
#####################################################



directory = '/mnt/gamma/user/sofya/scripts/final/0pipeline/13_LysEvo/codon_alignments/' 
files = os.listdir(directory)
# files = ['hom_cras_c3120_g1_i1.fasta']
for filename in files:
	codon_filename = filename[:-5] + 'nuc'
	tree_filename = filename[:-5] + 'nwk'
	tree_info, outdict = get_paml_tree(codon_filename, tree_filename)
	
	with open('./paml_ancestors/' + filename, 'w') as fout:
		fout.write('#' + tree_info + '\n')
		for key in outdict:
			fout.write('>' + key + '\n')
			fout.write(outdict[key] + '\n')


	# pause = int(input('procceed? ---> '))









# def parse_rst(path_to_rst):
	
# 	with open(path_to_rst, 'r') as fin:
# 		read_tree_info = False
# 		read_seqs = False
# 		for line in fin:
# 			if line == "tree with node labels for Rod Page's TreeView\n":
# 				read_tree_info = True 
# 				continue
# 			elif read_tree_info:
# 				tree_info = line.strip()
# 				read_tree_info = False
# 				continue
# 			elif line == "List of extant and reconstructed sequences\n":
# 				read_seqs = True
# 				continue
			
# 			elif read_seqs:
# 				if line[0] == '\n' or line[0] == ' ':
# 					continue

# 				if "Overall accuracy of the" in line:
# 					read_seqs = False
# 					break


# 				#################
# 				## U-hoo!    
# 				## look how many space in your rst does paml set for identificator (15 chars or more/less)
# 				## IF loop - use it if you don't want word "node" in your fasta
# 				#################
# 				if line[:4] == 'node':
# 					outdict[line[6:15].replace(" ", "")] = line[15:-1].replace(" ", "") # cuts off "node" - do you want it?
# 				else:
# 					outdict[line[:15].replace(" ", "")] = line[15:-1].replace(" ", "")



# 	return tree_info, outdict

