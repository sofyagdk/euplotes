#20-05-30 sofya 
# module dor alignment analyses 

#  s(n)- the weight of nucleotide N in a column, 
# si(n) = -ci(n) umnoj log(pi(n))
# ci - the amount of nucleotides N in a column i, 
# p -  the frequency of N in  column i. 
# Si = sum(by nucleotides n)si(n). 
# Important - counting p, it mustnt be 0 
# so add pseudocounts 1, for stability


import math 


def alignment_to_matrix(bp_dict):
	matrix = list()

	for species in bp_dict:
		matrix.append(list(bp_dict[species]))
	matrix = list(map(list, zip(*matrix))) #transposition - now column of the alignment corresponds to a line in the matrix
	return matrix

def count_column_ic(column):
	column_ic = 0 
	amount = {'A':0,  'T':0, 'C':0, 'G':0}
	exp = {'A':0.3,  'T':0.3, 'C':0.2, 'G':0.2}
	summ = 0

	for elem in column:
		if elem== '-':
			continue
		amount[elem] += 1 
		summ += 1
	if summ < 2:
		return 0 
	freq = dict()
	for nuc in amount:
		freq[nuc] = float(amount[nuc] + 0.25)/float(summ + 1)
		column_ic += freq[nuc]*math.log(freq[nuc]/exp[nuc], 2)
	
	
	return column_ic



def gap_cost(sequence, open_cost, add_cost):
	cost = 0
	in_indel = False
	for i in range(1, len(sequence)):
		if sequence[i] != '-':
			in_indel = False
			continue
		if in_indel== True:
			cost += add_cost
		else:
			in_indel== True
			cost += open_cost
	return cost 



def count_alignment_ic(bp_dict):
	alignment_ic = 0
	gaps = 0
	matrix = alignment_to_matrix(bp_dict)
	for column in matrix:
		alignment_ic += count_column_ic(column)
	for idr in bp_dict:
		gaps += gap_cost( bp_dict[idr], 0.2, 0.01)
	# print(alignment_ic)
	# print('GAPS', gaps)
	return alignment_ic - gaps 


def count_pairwise_score(seq1, seq2, match=1, mismatch = 1, gap_open = 1, gap_add = 0.1):
	score = 0 
	
	end = len(seq1) - 1
	while seq1[end] == '-':
		end  -= 1
	while seq2[end] == '-':
		end  -= 1

	i = 0
	while seq1[i] == '-':
		i += 1
	while seq2[i] == '-':
		i += 1
	
	in_indel = False

	while i<= end:
		if seq1[i] == seq2[i] and seq1[i] != '-' and seq2[i] != '-':
			score += match
		elif seq1[i] == '-' and seq2[i] == '-':
			i += 1
			continue
		elif seq1[i] != '-' and seq2[i] != '-' and seq1[i] == seq2[i]:
			score -= mismatch
		elif in_indel == False:
			score -= gap_open
			in_indel = True
		elif in_indel == True:
			score -= gap_add
		else:
			print('WHAT THE FUCK work betta batard')

		i+= 1
	# print(score)
	return score



def count_alignment_frompairwise(bp_dict):
	idrs = bp_dict.keys()
	alignment_score = 0 
	for i in range(len(idrs)):
		for j in range(i+1, len(idrs)):
			alignment_score += count_pairwise_score(bp_dict[idrs[i]], bp_dict[idrs[j]], 
				match=1, mismatch = 1, gap_open = 1, gap_add = 0.1)
	# print('sothefinal', alignment_score)
	return alignment_score