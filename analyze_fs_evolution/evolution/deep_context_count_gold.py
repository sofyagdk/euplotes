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


# filename = '/mnt/gamma/user/sofya/scripts/final/homo_bp_al_3/hom_cras_c24723_g1_i1.fasta'
directory ='/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/nuc_gold/'
# directory ='/mnt/gamma/user/sofya/scripts/final/0pipeline/nuc_alignments/'

files = os.listdir(directory)
# files = ['hom_foca_c28583_g5_i1.fasta']

# stranges = ['hom_foca_c13918_g2_i1.fasta', 'hom_foca_c30047_g1_i1.fasta',  'hom_foca_c31778_g1_i1.fasta', 'hom_eury_eury_c11198_g1_i1.fasta']

regexps  = dict() # key - the regexp line, first nom - is the amount when codon sencible, second nom - any codom orientation

# with open('./prev/deep_input_DEBUG.txt', 'r') as fin:
with open('./deep_input.txt', 'r') as fin:
	toalalyze = fin.readline()
	toalalyze = toalalyze.split()
	for elem in toalalyze:
		regexps[elem] = [0, 0]
		print(elem)


shis = 0
bef = 0 
AA = 0
AG = 0
conts_a = 0 
conts_g = 0
conts_r = 0

inserts_aa = 0 
inserts_ga = 0 
deletions_ata = 0 

# context = [CAGTRYSWKMBDHVN][CAGTRYSWKMBDHVN]A[AGR][ACM][ACM]

# contextmask = AAA[]AA


# AAAAAAAAAAAAAAAAAAch
# def makecodings(before_dict, species_list ):
# 	codings = dict()
# 	for idr in before_dict:
# 		coding, seq = mark.getcoding(boot.clean(before_dict[idr]), 35, 35, 20)
# 		if seq != boot.clean(before_dict[idr]):
# 			print('HOLY SHITTTT')

# 		if coding == list():
# 			return 'not_a_protein'
# 		new_coding = boot.get_inalignment_places(coding, before_dict[idr], boot.clean(before_dict[idr]))
		
# 		if idr[:4] in species_list:
# 			codings[idr[:4]] = new_coding
	
# 	return codings

def makecodings_alternative(before_dict, species_list ):
	codings = dict()
	for idr in before_dict:
		coding = fish_coding_alignment(idr)
		if idr[:4] in species_list:
			codings[idr[:4]] = coding
	
	return codings



def boarder(codings):
	left_boarders = list() 
	right_boarders = list()
	shifts = list()

	for sp in codings:
		left_boarders.append(codings[sp][0])
		right_boarders.append(codings[sp][-1])
		shifts.extend(codings[sp][1:-1])

	boarders = [max(left_boarders), min(right_boarders)]

	shiftset = set(shifts)
	shifts = list(shiftset)
	shifts.sort()
	return boarders, shifts

def cut_firstfew(seq):
	nom = [0, 1, 2]
	for i in range(0, len(seq), 3):
		if seq[i:i+3] in ['TAA', 'TAG']:
			nom.remove(0)
			break
	for i in range(1, len(seq), 3):
		if seq[i:i+3] in ['TAA', 'TAG']:
			nom.remove(1)
			break
	for i in range(2, len(seq), 3):
		if seq[i:i+3] in ['TAA', 'TAG']:
			nom.remove(2)
			break
	if len(nom) == 0:
		print(seq)
	return int(nom[0])



def get_search_region(seq, stencil_sequence,  stencil_sequence_idr, stencil_coding, boarders):
	search_region = str()
	coding = fish_coding_true(stencil_sequence_idr)



	shi = 2
	# print(coding )
	# print(stencil_coding)
	# print(stencil_sequence_idr)
	for i in range(1, len(stencil_coding) - 1, 2):
		shiftplace = stencil_coding[i]
		if shiftplace < boarders[1] and shiftplace > boarders[0]:
			if coding[shi] - coding[shi-1] == 1:
				stencil_sequence = stencil_sequence[:shiftplace] + '-' + stencil_sequence[shiftplace+1:]
			elif coding[shi] - coding[shi-1] == 2:
				stencil_sequence = stencil_sequence[:shiftplace] + '--' + stencil_sequence[shiftplace+2:]
		shi += 2

	for i in range(boarders[0], boarders[1]):
		if stencil_sequence[i] != '-':
			search_region += stencil_sequence[i]

	search_region = search_region[cut_firstfew(search_region[3:]):]

	return search_region


def count_contexts(search_region, context_regexp, context_len):
	
	context_expression = re.compile(context_regexp)
	amount = 0 

	for i in range(0, len(search_region), 3):
		string = search_region[i:i+context_len]
		if context_expression.match(string) != None:
			# print(string, context_expression, context_regexp)
			amount += 1
		# elif search_region[i+1:i+5] == 'AATA':
		# 	print('BUG!!!:   ', search_region[i+1:i+5], search_region[i-8:i+8] ,  context_regexp)
		# 	pause = int(input())

	return amount

def count_contexts_anyplace(search_region, context_regexp, context_len):
	
	context_expression = re.compile(context_regexp)
	amount = 0 

	for i in range(0, len(search_region), 1): # look here, not codon triples, but just 1!
		string = search_region[i:i+context_len]
		if context_expression.match(string) != None:
			amount += 1

	return amount


# def 
conts = 0 


# strangefile = open('strangefile_v8.txt', 'w')
# endsfile = open('endsfile_v8.txt', 'w')
# begsfile = open('begsfile_v8.txt', 'w')

# files = ['hom_cras_c11411_g1_i1.fasta']

for afile in files:
	# if afile in stranges:
	# 	continue

	print(afile)
	filename = directory + afile


	# before_dict = reader.readfasta_todictionary(filename)

	# bp_dict, matrix, species_list = gardener.turn_alignment_to_matrix(before_dict)
	# keys = list(before_dict.keys())
	# codings = makecodings(before_dict, species_list)
	# if codings == 'not_a_protein' or len(codings) == 0:
	# 	continue
	# ##########ALTERNATIVE##########
	# before_dict = reader.readfasta_todictionary(filename)

	# after_dict = dict()
	# for idr in before_dict:
	# 	coding, seq = mark.getcoding(boot.clean(before_dict[idr]), 35, 35, 20)
	# 	if seq != boot.clean(before_dict[idr]):
	# 		continue
	# 	else:
	# 		after_dict[idr] = before_dict[idr]

	# before_dict = dict(after_dict)

	# bp_dict, matrix, species_list = gardener.turn_alignment_to_matrix(before_dict)
	# keys = list(before_dict.keys())
	# codings = makecodings(before_dict, species_list)
	# if codings == 'not_a_protein' or len(codings) == 0:
	# 	continue

		# ##########ALTERNATIVE nom 2##########
	before_dict = reader.readfasta_todictionary(filename)
	idrs = before_dict.keys()

	bp_dict, matrix, species_list = gardener.turn_alignment_to_matrix(before_dict)
	for idr in idrs:
		if idr[0:4] ==species_list[0]:
			stencil_seq_idr = idr
	keys = list(before_dict.keys())
	codings = makecodings_alternative(before_dict, species_list)
	if codings == 'not_a_protein' or len(codings) == 0:
		continue


	# ###############################


	# ###############################

	boarders, shifts = boarder(codings)


	for sh in shifts:
		i = sh
		column = matrix[i]
		# local_tree = gardener.Ancestor_Tree('(petz:0.14339273,(raik:0.26772623,((octo:0.34300359,harp:0.21398805)12:0.08092917,(eury:0.21540617,(rari:0.27565106,(foca:0.26750949,(minu:0.12310897,cras:0.16946154)30:0.12093409)24:0.05101473)9:0.04135668)13:0.03839098)11:0.04916400):0.14339273);')
		# local_tree = gardener.Ancestor_Tree('(petz,(raik,((octo,harp),(eury,rari),(foca,(minu,cras)))));')
		local_tree = gardener.Ancestor_Tree('(petz:0.675090,(raik:1.138729,((octo:1.509930,harp:1.854564):0.176977,(eury:0.800680,(rari:1.121455,(foca:1.062187,(minu:0.721561,cras:0.730116):0.147409):0.155257):0.099063):0.084100):0.086484):0.409622);')
			
		local_tree.set_position(i)

		for j in range(len(column)):
			sp = species_list[j]
			coding = codings[sp]
			
			if sh in coding[1:-1]:
				m = 1
			else:
				m = 0 
			
			local_tree.get_leaves_by_name(sp)[0].node_mean(mean = m)
			
		meaningful_nodes = list()
		for n in local_tree.iter_leaves():
			if hasattr(n, 'mean') == True:
				meaningful_nodes.append(n.name)
		
		if len(meaningful_nodes) < 3:
			continue
		local_tree.prune(meaningful_nodes, preserve_branch_length = True)

		local_tree.node_mean()
		local_tree.node_mean_downgoing()
		# local_tree.view_means()

		#useful stuff here! no deleting
		if local_tree.mean == 1:
			shis += 2
		
		for n in local_tree.iter_descendants():
			if n.is_leaf() == True:
				continue
			if n.mean == 1:
				shis += 2
		print(shis)





	i = (boarders[1] + boarders[0])//2
	c = (boarders[1] - boarders[0])//2

	column = matrix[i]
	# local_tree = gardener.Ancestor_Tree('(petz:0.14339273,(raik:0.26772623,((octo:0.34300359,harp:0.21398805)12:0.08092917,(eury:0.21540617,(rari:0.27565106,(foca:0.26750949,(minu:0.12310897,cras:0.16946154)30:0.12093409)24:0.05101473)9:0.04135668)13:0.03839098)11:0.04916400):0.14339273);') 
	# local_tree = gardener.Ancestor_Tree('(petz,(raik,((octo,harp),(eury,rari),(foca,(minu,cras)))));')
	local_tree = gardener.Ancestor_Tree('(petz:0.675090,(raik:1.138729,((octo:1.509930,harp:1.854564):0.176977,(eury:0.800680,(rari:1.121455,(foca:1.062187,(minu:0.721561,cras:0.730116):0.147409):0.155257):0.099063):0.084100):0.086484):0.409622);')
			
	local_tree.set_position(i)

	stencil_coding = codings[species_list[0]]
	stencil_seq = bp_dict[species_list[0]]
	# print('the selectedd animal', species_list[0])


	for j in range(len(column)):
		sp = species_list[j]
		context_string = gardener.cut_context(i, sp, c, bp_dict)
			
		local_tree.get_leaves_by_name(sp)[0].leaf_context(context_string)
			
	meaningful_nodes = list()
	for n in local_tree.iter_leaves():
		if hasattr(n, 'context') == True:
			meaningful_nodes.append(n.name)
		
	if len(meaningful_nodes) < 3:
		continue
	
	local_tree.prune(meaningful_nodes, preserve_branch_length = True)
	local_tree.node_context_upgoing(2*c)
	local_tree.node_context_downgoing(2*c)



	# local_tree.view_contexts()
	#conts += count_contexts(string_tosearch, stencil_sequence, regexp)
	print(len(local_tree.context), local_tree.context)
	search_region = get_search_region(local_tree.context, stencil_seq, stencil_seq_idr, stencil_coding, boarders)
	print(len(search_region ))
	# pause = int(input())


	for key in regexps:
		regexps[key][0] += 2*count_contexts(search_region, key, 6)
		regexps[key][1] += 2*count_contexts_anyplace(search_region, key, 6)

	# conts_a += 2*count_contexts(search_region, 'AAAA[GA]', 5)
	# conts_g += 2*count_contexts(search_region, 'AAGA[AG]', 5)
	# # conts_r += count_contexts(search_region, '[CAGTRYSWKMBDHVN][CAGTRYSWKMBDHVN]A[R][ACM][ACM]', 6)
	# # conts_a += count_contexts(search_region, '[CAGTRYSWKMBDHVN][CAGTRYSWKMBDHVN]A[R][ACM][ACM]', 6)
	# inserts_aa += 2*count_contexts(search_region, '[CAGTRYSWKMBDHVN][CAGTRYSWKMBDHVN]AA', 4) 
	# inserts_ga += 2*count_contexts(search_region, '[CAGTRYSWKMBDHVN][CAGTRYSWKMBDHVN]A[G]', 4) 
	# deletions_ata += 2*count_contexts(search_region, '[CAGTRYSWKMBDHVN][CAGTRYSWKMBDHVN]A[T][A]', 5) 


	# print(local_tree.position)
	# print(boarders, shifts)
	# print(local_tree.context)
	# print(local_tree)
	# print(afile)


	# print('----------------------------------------------------')

print('a', conts_a)
print('g', conts_g)
print('r', conts_r)
print(shis)


print('insaa', inserts_aa) 
print('insga', inserts_ga) 
print('del ata', deletions_ata)

with open('./contexts_forced_treehairbrush_gold.txt', 'a') as fout:
	for key in regexps:
		out = [key, str(regexps[key][0]), str(regexps[key][1]), '\n']
		fout.write('\t'.join(out))

