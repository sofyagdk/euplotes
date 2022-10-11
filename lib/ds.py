import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader
import math
translation_table = {
"TTT" : "F", "TCT" : "S", "TAT" : "Y", "TGT" : "C", 
"TTC" : "F", "TCC" : "S", "TAC" : "Y", "TGC" : "C", 
"TTA" : "L", "TCA" : "S", "TAA" : "*", "TGA" : "C", 
"TTG" : "L", "TCG" : "S", "TAG" : "*", "TGG" : "W", 
"CTT" : "L", "CCT" : "P", "CAT" : "H", "CGT" : "R", 
"CTC" : "L", "CCC" : "P", "CAC" : "H", "CGC" : "R", 
"CTA" : "L", "CCA" : "P", "CAA" : "Q", "CGA" : "R", 
"CTG" : "L", "CCG" : "P", "CAG" : "Q", "CGG" : "R", 
"ATT" : "I", "ACT" : "T", "AAT" : "N", "AGT" : "S", 
"ATC" : "I", "ACC" : "T", "AAC" : "N", "AGC" : "S", 
"ATA" : "I", "ACA" : "T", "AAA" : "K", "AGA" : "R", 
"ATG" : "M", "ACG" : "T", "AAG" : "K", "AGG" : "R", 
"GTT" : "V", "GCT" : "A", "GAT" : "D", "GGT" : "G", 
"GTC" : "V", "GCC" : "A", "GAC" : "D", "GGC" : "G", 
"GTA" : "V", "GCA" : "A", "GAA" : "E", "GGA" : "G", 
"GTG" : "V", "GCG" : "A", "GAG" : "E", "GGG" : "G" 
}



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


aminoacids = {'Y': ['TAT', 'TAC'], 'C': ['TGC', 'TGT', 'TGA'], 
'F': ['TTT', 'TTC'], 'I': ['ATC', 'ATT', 'ATA'], 
'L': ['CTC', 'CTT', 'CTG', 'CTA', 'TTG', 'TTA'], 
'K': ['AAG', 'AAA'], 'H': ['CAC', 'CAT'], 'G': ['GGT', 'GGC', 'GGG', 'GGA'], 
'D': ['GAC', 'GAT'], 'V': ['GTT', 'GTG', 'GTA', 'GTC'], 'X': ['TAG', 'TAA'], 
'P': ['CCC', 'CCG', 'CCA', 'CCT'], 'A': ['GCC', 'GCG', 'GCA', 'GCT'], 
'W': ['TGG'], 'M': ['ATG'], 'Q': ['CAA', 'CAG'],
 'R': ['AGA', 'CGG', 'CGC', 'CGT', 'AGG', 'CGA'], 
 'S': ['TCC', 'TCT', 'TCG', 'AGT', 'AGC', 'TCA'], 'N': ['AAC', 'AAT'], 
 'E': ['GAA', 'GAG'], 'T': ['ACT', 'ACC', 'ACG', 'ACA']}

d4 = ['GGT', 'GGC', 'GGG', 'GGA',
 'GTT', 'GTG','GTA', 'GTC', 
 'CCC', 'CCG', 'CCA', 'CCT', 
 'GCC', 'GCG', 'GCA', 'GCT', 
 'ACT', 'ACC', 'ACG', 'ACA']


bp = ['A', 'T', 'G', 'C']

#1 we need to create a function that reads a string an for every codon writes the synonymous 

def different(b):
	a = ['A', 'T', 'G', 'C']
	a.remove(b)
	return a


def synonym(codon):
	aa = translation_table[codon]
	s = 0 
	ch = list()
	
	for i in [0, 2]:
		for b in different(codon[i]):
			ch.append(list(codon))
			ch[-1][i] = b
	for change in ch:
		newa = translation_table[''.join(change)]
		if newa == aa:
			s += 1

	return float(s)/float(3)  




def issyn(c1, c2):
	if translation_table[''.join(c1)] == translation_table[''.join(c2)]:
		return True
	else:
		return False


def substitutions(c1, c2):
	dis = 0
	s = 0 
	diffs = list()
	
	for i in range(3):
		if c1[i] != c2[i]:
			dis += 1
			diffs.append(i)

	if dis == 0:
		return 0, 0


	if dis == 1:
		return 1, int(issyn(c1, c2))

	if dis == 2:
		for d in diffs:
			z1 = list(c1)
			z1[d] = c2[d]
			s += int(issyn(c1, z1)) + int(issyn(z1, c2))
		return 2, s/2.0

	if dis == 3:
		for d1 in range(3):
			z1 = list(c1)
			z1[d1] = c2[d1]
			for d2 in range(3):
				if d1 == d2:
					continue
				z2 = list(z1)
				z2[d2] = c2[d2]
				s += int(issyn(c1, z1)) + int(issyn(z1, z2)) + int(issyn(z2, c2))
		
		return 3, s/6.0

# print(substitutions(list('CAC'), list('TGG')))
# print(substitutions(list('CAC'), list('TGG')))

# exit()
syn = dict()
non = dict()

for c1 in codon_list:
	syn[c1] =  dict()
	non[c1] =  dict()
	for c2 in codon_list:
		dist, s = substitutions(list(c1), list(c2))
		# print(s, end = ' ')
		syn[c1][c2] = s
		non[c1][c2] = dist - s



def mydnds(seq1, seq2):


	s_sum = 0 
	

	sd = 0 
	nd = 0 

	for i in range(0, len(seq1), 3):
		c1 = seq1[i:i+3]
		c2 = seq2[i:i+3]
		
		s_sum += synonym(c1)
		s_sum += synonym(c2)

		sd += syn[c1][c2]
		nd += non[c1][c2]


	r = len(seq1)//3
	S = float(s_sum)/2.0
	N = 3*r - S


	return nd, N, sd, S


def diff(seq1, seq2):
	d = 0 
	for i in range(len(seq1)):
		if seq1[i] != seq2[i]:
			d += 1 
	return d 



def cleanfromgaps(align_list):
	gap_l = 0 
	for e in align_list[0]:
		if e == '-':
			gap_l += 1 
		else:
			if gap_l%3 != 0:
				print("check/ds.py:: cleanfromgaps() function:: ERROR GAP NOT LENGTH 3")
				pause = int(input())
	for e in align_list[1]:
		if e == '-':
			gap_l += 1 
		else:
			if gap_l%3 != 0:
				print("check/ds.py:: cleanfromgaps() function:: ERROR GAP NOT LENGTH 3")
				pause = int(input())



	i = 0 
	while i < len(align_list[0]):
		if align_list[0][i] == '-':
			align_list[0] = align_list[0][:i] + align_list[0][i + 1:]
			align_list[1] = align_list[1][:i] + align_list[1][i + 1:]
		elif align_list[1][i] == '-':
			align_list[0] = align_list[0][:i] + align_list[0][i + 1:]
			align_list[1] = align_list[1][:i] + align_list[1][i + 1:]
		else:
			i += 1
	return align_list


def getdnds(seq1, seq2):
	s1, s2 = cleanfromgaps([seq1, seq2])
	nd, N, sd, S  = mydnds(s1, s2)
	return nd, N, sd, S 


