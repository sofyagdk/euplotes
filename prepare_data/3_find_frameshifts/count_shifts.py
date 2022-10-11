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

with open('/mnt/gamma/user/sofya/scripts/final/0pipeline/18_graphmark/mark_raw.txt', 'r') as fin:
	for line in fin:
		cod = line.split()[2]
		cod = list(map(int, cod.split(',')))
		codings_true[line.split()[1]] = cod


# codings_alignment = dict()

# with open('/mnt/gamma/user/sofya/scripts/final/0pipeline/18_graphmark/mark_raw.txt', 'r') as fin:
# 	for line in fin:
# 		cod = line.split()[4]
# 		cod = list(map(int, cod.split(',')))
# 		codings_true[line.split()[1]] = cod


###########################################
### Create dictionaries
###########################################


for idr in codings_true:
	fs_1 = 0
	fs_2 = 0
	rthr = 0
	coding  = codings_true[idr]
	for i in range(2, len(coding), 2):
		s = coding[i]
		if coding[i] - coding[i-1] == 1:
			fs_1 += 1
		elif coding[i] - coding[i-1] == 2:
			fs_2 += 1
	print(' '.join([idr, str(fs_1), str(fs_2)]))
