import os
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader
import evo_ans_temp_v4 as gardener
from ete3 import Tree





#####################################################
### Dictionaries and Data
#####################################################


codings_true = dict()

with open('/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/mark.txt', 'r') as fin:
	for line in fin:
		cod = line.split()[2]
		cod = list(map(int, cod.split(',')))
		codings_true[line.split()[1]] = cod

newickFullTree = '(petz:0.14339273,(raik:0.26772623,\
	((octo:0.34300359,harp:0.21398805)12:0.08092917,\
	(eury:0.21540617,(rari:0.27565106,(foca:0.26750949,\
	(minu:0.12310897,cras:0.16946154)30:0.12093409)\
	24:0.05101473)9:0.04135668)13:0.03839098)11:0.04916400)R\
	:0.14339273);'

	# dS tree:
newickFullTree = '(petz: 0.675090, (raik: 1.138729, ((octo: 1.509930, harp: 1.854564): 0.176977, (eury: 0.800680, (rari: 1.121455, (foca: 1.062187, (minu: 0.721561, cras: 0.730116): 0.147409): 0.155257): 0.099063): 0.084100): 0.086484): 0.409622);'




#####################################################
### Functions
#####################################################


def write_nuc(alignment_dict, filename):
	with open(filename, 'w') as fout:
		fout.write('\t' + str(len(alignment_dict))+ '\t' + str(len(alignment_dict[alignment_dict.keys()[0]])) + '\n')
		for key in alignment_dict:
			animal = key[:4]
			fout.write(animal + '\n')
			fout.write(alignment_dict[key] + '\n')

def write_tree(alignment_dict, filename):
	local_tree = gardener.Ancestor_Tree( newickFullTree, format = 1)
	leaves = list()
	for idr in alignment_dict.keys():
		animal = idr[:4]
		leaves.append(animal)
	local_tree.prune(leaves, preserve_branch_length = True)
	with open(filename, 'w') as fout:
		fout.write('\t' + str(len(leaves))+ '\t1\n')
		# format = 5 means saving branch lengths and NOT displaying internal node names
		fout.write(local_tree.write(format = 5)) 


#####################################################
### DATA
#####################################################
directory = '/mnt/gamma/user/sofya/scripts/final/0pipeline/13_LysEvo/codon_alignments/' 
files = os.listdir(directory)




#####################################################
### Main part
#####################################################

already = os.listdir('/mnt/gamma/user/sofya/scripts/final/0pipeline/13_LysEvo/codon_alignments_nuc/')

for filename in files:
	if filename[:-5] + 'nuc' in already:
		continue

	path = directory + filename
	bp_dict = reader.readfasta_todictionary(path)
	nucfilename = '/mnt/gamma/user/sofya/scripts/final/0pipeline/13_LysEvo/codon_alignments_nuc/' + filename[:-5] + 'nuc'
	write_nuc(bp_dict, nucfilename)
	nwkfilename = '/mnt/gamma/user/sofya/scripts/final/0pipeline/13_LysEvo/nwk_trees/'+ filename[:-5] + 'nwk'
	write_tree(bp_dict, nwkfilename)




