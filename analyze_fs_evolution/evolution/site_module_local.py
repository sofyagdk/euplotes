#november 17 with new finder
import sys

sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import evo_ans_temp_v4 as gardener
from ete3 import Tree
import reader
import os 
import subprocess
# import marking_module_v3 as mark
import bootstrap_tables_v4 as boot 
# import marker_pnpsbased as something_new

# filename = '/mnt/gamma/user/sofya/scripts/final/homo_bp_al_3/hom_cras_c24723_g1_i1.fasta'
directory ='/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/nuc_gold/'
files = os.listdir(directory)


# del_file = open('deletions_0808.txt', 'w')
# ins_file = open('insertions_0808.txt', 'w')

codingsfilename = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/mark.txt'

# self.string_branch(self.make_branch_that_leads_to(n))
def fish_coding(idr):
	global codingsfilename
	with open(codingsfilename, 'r') as fin:
		for line in fin:
			linelist = line.split()
			# cras_c3120_g1_i1	[35, 842, 843, 885]	[40, 848, 849, 920]	[40, 848]
			if linelist[1] == idr:
				cod = linelist[4]
				cod = list(map(int, cod.split(',')))
				return cod



def makecodings(before_dict, species_list ):
	codings = dict()
	k = before_dict.keys()
	ref = k[0]

	for idr in before_dict:
		# new_coding = boot.get_inalignment_places(coding, before_dict[idr], boot.clean(before_dict[idr]))
		new_coding = fish_coding(idr)
		if idr[:4] in species_list:
			codings[idr[:4]] = new_coding
	return codings

# def boarder(codings):
# 	left_boarders = list() 
# 	right_boarders = list()
# 	shifts = list()

# 	for sp in codings:
# 		left_boarders.append(codings[sp][0])
# 		right_boarders.append(codings[sp][-1])
# 		shifts.extend(codings[sp][1:-1])

# 	boarders = [max(left_boarders), min(right_boarders)]
# 	# print(shifts)
# 	shiftset = set(shifts)
# 	shifts = list(shiftset)
# 	shifts.sort()
# 	return boarders, shifts

def boarder(codings):
	left_boarders = list() 
	right_boarders = list()
	shifts = list()

	print(codings)
	for sp in codings:
		left_boarders.append(codings[sp][0])
		right_boarders.append(codings[sp][-1])
		shifts.extend(codings[sp][1:-1])

	boarders = [min(left_boarders), max(right_boarders)]
	# print(shifts)
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
	
	return int(nom[0])

# def get_search_region(seq, stencil_sequence, stencil_coding, boarders):
# 	search_region = str()
# 	coding, seq_sten = mark.getcoding(boot.clean(stencil_sequence), 35, 35, 20)
# 	coding = something_new.get_coding_pnps_consult(idr, k[1])

# 	shi = 2
# 	for shiftplace in stencil_coding[1:-1]:
# 		if shiftplace < boarders[1] and shiftplace > boarders[0]:
# 			if coding[shi] - coding[shi-1] == 1:
# 				stencil_sequence = stencil_sequence[:shiftplace] + '-' + stencil_sequence[shiftplace+1:]
# 			elif coding[shi] - coding[shi-1] == 2:
# 				stencil_sequence = stencil_sequence[:shiftplace] + '--' + stencil_sequence[shiftplace+2:]
# 		shi += 2

# 	for i in range(boarders[0], boarders[1]):
# 		if stencil_sequence[i] != '-':
# 			search_region += stencil_sequence[i]

# 	search_region = search_region[cut_firstfew(search_region):]

# 	return search_region

def codon_place(seq, coding, position):
	last = 0 
	while coding[last] <= position:
		last += 1
	last -= 1

	cut = boot.clean(seq[coding[last]:position])

	return len(cut)%3



context_file = open('events_treehairbrush_golden.txt', 'w')
# context_file = open('events_tree2_golden.txt', 'w')

context_file.write('event_type	ancestor_context	child_context	place	alignment_position	codon_position	start	end	procent	branch	filename\n')

# files = ['hom_cras_c4600_g1_i1.fasta']


for filename in files:
	print(filename)
	path = directory + filename
	before_dict = reader.readfasta_todictionary(path)
	before_dict = reader.filter_paralogues(before_dict)


	bp_dict, matrix, species_list = gardener.turn_alignment_to_matrix(before_dict)
	keys = list(before_dict.keys())
	codings = makecodings(before_dict, species_list)
	# print(codings)
	# pause = int(input())
	if codings == 'not_a_protein' or len(codings) == 0:
		continue

	boarders, shifts = boarder(codings)

	stencil_coding = codings[species_list[0]]
	stencil_seq = bp_dict[species_list[0]]

	for i in range(5,len(matrix)-5):
		
		column = matrix[i]
		if '-' in column:
			if len(set(column)) == 1:
				continue
			# local_tree = gardener.Ancestor_Tree('(petz:0.14339273,(raik:0.26772623,((octo:0.34300359,harp:0.21398805)12:0.08092917,(eury:0.21540617,(rari:0.27565106,(foca:0.26750949,(minu:0.12310897,cras:0.16946154)30:0.12093409)24:0.05101473)9:0.04135668)13:0.03839098)11:0.04916400):0.14339273);')
			# local_tree = gardener.Ancestor_Tree('(petz,(raik,((octo,harp),(eury,rari),(foca,(minu,cras)))));')
			local_tree = gardener.Ancestor_Tree('(petz:0.675090,(raik:1.138729,((octo:1.509930,harp:1.854564):0.176977,(eury:0.800680,(rari:1.121455,(foca:1.062187,(minu:0.721561,cras:0.730116):0.147409):0.155257):0.099063):0.084100):0.086484):0.409622);')
			
		else:
			continue
		# onebpgap = False
		local_tree.set_position(i)

		for j in range(len(column)):
			sp = species_list[j]
			context_string = gardener.cut_context(i, sp, 5, bp_dict)
			
			if column[j] == '-':
				# if context corresponds to somethig - m = 0, otherwise continue
				# print('jskhfahflkajhf', column)
				# if context_string[2] != '-' and context_string[4] != '-':
				# 	onebpgap = True
				m = 0
			# elif column[j] == 'T':
			# 	m = 1 
			else:
				m = 1
				# continue

			
			# if context_string[len(context_string)//2] == '-' and context_string[len(context_string)//2 - 1] != '-' and context_string[len(context_string)//2+1] != '-':
			if context_string[len(context_string)//2 - 1:] != '-' and context_string[len(context_string)//2+1] != '-':
				local_tree.get_leaves_by_name(sp)[0].node_mean(mean = m)
				local_tree.get_leaves_by_name(sp)[0].leaf_context(context_string)
		meaningful_nodes = list()
		for n in local_tree.iter_leaves():
			if hasattr(n, 'mean') == True:
				meaningful_nodes.append(n.name)
		if len(meaningful_nodes) < 3:
			print('too little')
			continue
		
		local_tree.prune(meaningful_nodes, preserve_branch_length = True)
		
		if local_tree.is_onebp_gap() == False:
			with open('mistake_1.txt', 'a') as fout:
				fout.write(filename + '\t' +  str(i) + '\n')
				# pause = input()
			continue


		if local_tree.position < boarders[0] or local_tree.position > boarders[1]:
			place = 'out'
			incodon = -1
		elif local_tree.position in shifts:
			place = 'shi'
			# incodon = codon_place(stencil_seq, stencil_coding, local_tree.position)
		elif local_tree.position == boarders[0] or local_tree.position == boarders[1]:
			place = 'bor'	
			# incodon = codon_place(stencil_seq, stencil_coding, local_tree.position)
		else:
			place = 'ins'
			# incodon = codon_place(stencil_seq, stencil_coding, local_tree.position)
		

		local_tree.node_mean()
		local_tree.node_mean_downgoing()
		local_tree.view_means()

		local_tree.node_context_upgoing(11)
		local_tree.node_context_downgoing(11)
		local_tree.view_contexts()


		local_tree.event_nodes_maker(local_tree)
		event_nodes = local_tree.event_nodes

		print(event_nodes)

		print(local_tree.position)
		print(filename)

		try:
			procent = float(len(boot.clean(stencil_seq[boarders[0]:local_tree.position])))/float(len(boot.clean(stencil_seq[boarders[0]:boarders[1]])))
		except ZeroDivisionError:
			procent = 0

		for node in event_nodes['deletion']:
			ancestor_context = node._up.context
			child_context = node.context
			# outline = ['deletion', ancestor_context, child_context, place, local_tree.position, incodon, boarders[0], boarders[1], procent, local_tree.string_branch(local_tree.make_branch_that_leads_to(node)), afile, '\n']
			outline = ['deletion', ancestor_context, child_context, place, local_tree.position, 0, boarders[0], boarders[1], procent, local_tree.string_branch(local_tree.make_branch_that_leads_to(node)), filename, '\n']
			context_file.write('\t'.join(list(map(str, outline))))
		for node in event_nodes['insertion']:
			ancestor_context = node._up.context
			child_context = node.context
			# outline = ['insertion', ancestor_context, child_context, place, local_tree.position, incodon, boarders[0], boarders[1], procent, local_tree.string_branch(local_tree.make_branch_that_leads_to(node)), afile, '\n']
			outline = ['insertion', ancestor_context, child_context, place, local_tree.position, 0, boarders[0], boarders[1], procent, local_tree.string_branch(local_tree.make_branch_that_leads_to(node)), filename, '\n']
			context_file.write('\t'.join(list(map(str, outline))))
		

		print('----------------------------------------------------')



context_file.close()
