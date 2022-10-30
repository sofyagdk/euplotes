from __future__ import print_function
import os
import sys
sys.path.append("./../../../../lib/")
import reader
import ds
import subprocess
import glob
import random
import math 
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def get_markup(idr):
	with open('./../../mark.txt', 'r') as fin:
		for line in fin:
			if idr in line:
				copy = line.split()
				if idr == copy[1]:
					return line[:-1]

# 3	cras_c6209_g1_i1	100.000	30	0	0	1	30	558	587	5.66e-10	56.5

def get_ribobits(idr):
	with open('mapcrassus.txt', 'r') as finmap:
		for line in finmap:
			copy = line.split()
			if copy[0] != idr:
				continue
			if int(copy[2]) != 0:
				return list()
			idr, straight, reverse, places = copy
			return places
	return list()

def modify_ribobits(ribobits):
	ribobits = ribobits.split(',')
	ribobits = ribobits[:-1] 
	newbits = list()
	for i in range(0, len(ribobits), 2):
		s = int(ribobits[i])
		e = int(ribobits[i+1])
		newbits.append([s, e])
	newbits.sort(key=lambda x:x[0])
	return newbits

# CTCTCAAGAAGAGAGAACTACCTTAAAATT
# 3	cras_c6209_g1_i1	100.000	30	0	0	1	30	558	587	5.66e-10	56.5

directory = '../../../data/nuc_alignments/'
files = os.listdir(directory)


for filename in files:

	path = directory + filename
	bp_dict = reader.readfasta_todictionary(path)
	bp_dict = reader.filter_paralogues(bp_dict)
	idrs = bp_dict.keys()

	for idr in bp_dict:
		if 'cras' in idr:
			seq = reader.clean(bp_dict[idr])
			ribobits = get_ribobits(idr)
			if ribobits == list():
				continue
			ribobits = modify_ribobits(ribobits)
			newseq = [0]*len(seq)
			for ribo_tuple in ribobits:
				s = ribo_tuple[0]
				e = ribo_tuple[1]
				for t in range(s-1, e): # see CTCTCAAGAAGAGAGAACTACCTTAAAATT in blast:3	cras_c6209_g1_i1	100.000	30	0	0	1	30	558	587	5.66e-10	56.5
					newseq[t] += 1
			info = get_markup(idr)
			print('>' + idr + ' ' + info)
			print(' '.join(list(map(str, newseq))))
# 			cut, has_shift = getcut(bp_dict, idr, filename)
