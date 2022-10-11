# AATA
# AGTA
# # 


import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
sys.path.append("/mnt/gamma/user/sofya/scripts/final/0pipeline/13_LysEvo")
import evo_ans_temp_v4 as gardener
from ete3 import Tree
import reader
import os 
import subprocess
# import marking_module_v3 as mark
# import bootstrap_tables_v4 as boot 
import re
from string_tree_class import String_Tree

###########################################
### Create dictionaries
###########################################

codings_true = dict()

with open('/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/mark.txt', 'r') as fin:
	for line in fin:
		cod = line.split()[2]
		cod = list(map(int, cod.split(',')))
		codings_true[line.split()[1]] = cod


# codings_alignment = dict()

# with open('/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/mark.txt', 'r') as fin:
# 	for line in fin:
# 		cod = line.split()[4]
# 		cod = list(map(int, cod.split(',')))
# 		codings_true[line.split()[1]] = cod


###########################################
### Functions
###########################################


def clean_paml_tree_line(tree_line):
	# tree_line = tree_line.replace('_', '')
	splited = tree_line.split('_')
	for i in range(len(splited) - 1):
		splited[i] = splited[i][:-1]
	return ''.join(splited)


def branched_and_named(tree_named, tree_branched):
	tree_new = str()


	TN = tree_named.split()
	nodenames = list()
	for n in range(len(TN)):
		if TN[n][-1] == ')':
			nodenames.append(TN[n+1])

	NOW = 0 
	tree_new = str()
	b = 0 
	while b < len(tree_branched):
		if tree_branched[b:b+3] == ')1:':
			tree_new += ')' + nodenames[NOW] + ':'
			NOW += 1
			b += 3
		elif tree_branched[b] == ';':
			tree_new += nodenames[NOW] + ';'
			b += 1

		else:
			tree_new += tree_branched[b]
			b += 1

	return tree_new


def norm_alignment(alt_seq, position):
	filled = -1 
	now = -1 
	while  (filled < position):
		now += 1
		if alt_seq[now] != '-':
			filled += 1
	return now



# def cuttoboarder(bp_dict):
# 	beg = 0
# 	end = len(bp_dict[bp_dict.keys()[0]])
# 	for elem in bp_dict:
# 		seq = bp_dict[elem]

# 		for i in range(len(seq)):
# 			if seq[i] != '-':
# 				beg = max(beg, i)
# 				break

# 		for i in range(1, len(seq)):
# 			if seq[-i] != '-':
# 				end = min(end, len(seq) - i + 1)
# 				break
# 	beg = beg + (3 - beg%3)
# 	end = end + (3 - end%3)

# 	print(beg, end, len(bp_dict[bp_dict.keys()[0]]))

# 	return beg, end






###########################################
### Main Part 
###########################################




# count the branches length 

# 1) AAA_AR in parent node 
# 2) AAR_AR in parent node 
# 3) AR in parent node 
# 4) anything in parent node

# all must not be the ones that stand in the listdir

# function that takes another function - checks the node to be right 

# tree with parental sequences in nodes!

# any type of passage and getting the branch lengths





BUG = 0


paml_reconstruction_directory = '/mnt/gamma/user/sofya/scripts/final/0pipeline/13_LysEvo/paml_ancestors/'
nucleotide_directory = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/nuc_gold/'
files = os.listdir(paml_reconstruction_directory)




for filename in files:
	# if filename not in [ 'hom_cras_c26056_g1_i1.fasta' ]:
	# 	continue
	# pause = int(input())

	#this is because im a lazy idiot but who cares
	path = nucleotide_directory + filename
	bp_dict = reader.readfasta_todictionary(path)
	identificators = bp_dict.keys()



	# actually main part

	###################################
	###### Tree reading and initiation
	###################################

	path = paml_reconstruction_directory + filename
	bp_dict = reader.readfasta_todictionary(path)
	
	if len(bp_dict) == 3:
		continue
	with open(path, 'r') as fin:
		infoline = fin.readline()[1:-1]
	tree_line = clean_paml_tree_line(infoline)



	meaningful_nodes = list()
	for idr in identificators:
		meaningful_nodes.append(idr[:4])
	# local_tree = gardener.Ancestor_Tree('(petz: 0.180584, (raik: 0.180584, ((octo: 0.117206, harp: 0.117206): 0.046075, (eury: 0.147157, (rari: 0.126699, (foca: 0.103248, (minu: 0.066474, cras: 0.066474): 0.036774): 0.023451): 0.020458): 0.016124): 0.017303): 0.000000);',
	# format = 1)
	local_tree = gardener.Ancestor_Tree('(petz:0.675090,(raik:1.138729,((octo:1.509930,harp:1.854564):0.176977,(eury:0.800680,(rari:1.121455,(foca:1.062187,(minu:0.721561,cras:0.730116):0.147409):0.155257):0.099063):0.084100):0.086484):0.409622);',
		format = 1)
	
	local_tree.prune(meaningful_nodes, preserve_branch_length = True)
	tree_size = local_tree.tree_distance()

	tree_named = tree_line
	tree_branched = local_tree.write()
	tree_new = branched_and_named(tree_named, tree_branched)


	


	shift_positions = list()
	# boarders = list()



	local_tree = String_Tree(tree_new, format = 1)

	# local_tree.pseudo_init()

	for elem in bp_dict:
		shift_list = list()
		nodename = elem       
		if elem[:4] == 'node':
			nodename = elem[5:]
		else:
			
			cod = 0
			already_passed = 0


			for idr in identificators:
				if idr[:4] == elem:
					cod = codings_true[idr]
					# boarders.append([norm_alignment(bp_dict[elem],cod[0]),  norm_alignment(bp_dict[elem],cod[-1])])
					break
			
			
			for i in range(2, len(cod), 2):

				shift_list.append( norm_alignment(bp_dict[elem], already_passed + cod[i-1] - cod[i-2]))
				# print(shift_list)
				already_passed  += cod[i-1] - cod[i-2]

		shift_positions.extend(shift_list)


		node = local_tree.search_nodes(name=nodename)[0]
		node.set_sequence(bp_dict[elem])
		node.set_shifts(shift_list)



	shift_positions = list(set(shift_positions))
	local_tree.solve_shifts()





 	# ancestor_context 	child_context 	place 	alignment_position 	codon_position 	start 	end 	procent 	branch 	filename 	inzoo
	# 28 	deletion 	GTAAATAACAT 	GTAAA-AACAT 	shi 	452 	0 	38 	1035 	0.417178 	rari,eury,petz,foca,octo_minu,cras 	hom_cras_c29146_g1_i1.fasta 	True
	# 162 	deletion 	TGAAATAAAAC 	TGAAA-AAGAC 	shi 	157 	0 	1 	1398 	0.104820 	octo,raik,petz_eury,rari 	hom_eury_eury_c19346_g1_i1.fasta 	True
	# 316 	deletion 	ATAAATAAAAA 	ATAAA-AGGAA 	shi 	220 	0 	1 	637 	0.331190 	rari,octo,minu,foca_cras 	hom_cras_c5844_g1_i1.fasta 	True
	# 407 	deletion 	TTAAATAAGAC 	TGAAA-AAGAC 	shi 	300 	0 	1 	755 	0.419308 	eury,foca,cras_rari 	hom_cras_c4123_g1_i1.fasta 	True




	# tree_size = local_tree.tree_distance()
	ancestor_sequence = local_tree.sequence
	ancestor_shifts = local_tree.shift_events


	# cheking for no "-" in +- 5 
	# registaring changes filled -> gapped = loss, gapped -> filled = gain 
	# and 
	# before + insertion + after 
	# insertion_codon_place = [0, 1, 2]
	# filled = before + insertion + after 
	# gapped = before + ('-' * len(insertion)) + after 
	
	codon_position = 1
	tocheck = 'AATA'
	tocheck = 'AAA'

	branchlen = 0
	number = 0
	# notinserted = 'AAA' #-> Naa_a


	inRoot = True
	# beg, end = cuttoboarder(bp_dict)


	for node in local_tree.traverse("preorder"): 
		# print(node)
		# print(node.shift_events)
		if node.is_leaf():
			continue
		if inRoot:
			inRoot = False
			continue
		if node.is_root():
			continue







		seq_len = len(node.sequence)

		for i in range(0, seq_len-3, 3):

			seq_ancestor = node.sequence[i + codon_position:i+codon_position+len(tocheck)]



			if i+3 in node.shift_events:
				continue


			for child in node._children:

				seq_child = child.sequence[i + codon_position:i+codon_position+len(tocheck)]
				child_dist = child.get_distance(node)




				################ strict alternative ###############

				
				if ((seq_child == tocheck) and (seq_ancestor == tocheck)):
					branchlen += child_dist
					number += 1
					# BUG += number
					# print(seq_ancestor, seq_child, child_dist, i%3)




				# ############ non-strict alternative #############

				# if seq_ancestor == tocheck:
				# 	branchlen += child_dist
				# 	print(seq_ancestor)

				# ################# end alternative #################





	print(' '.join(list(map(str, [filename, seq_len, tocheck, codon_position,  branchlen, number ]))))



	# loss = 0
	# gain = 0
	# B_shift = 0 
	# B_AAA = 0 
	# B_AAR = 0
	# B_AR = 0
	# B_anything = 0 

# print(BUG)