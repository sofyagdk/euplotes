

import os
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader
import random
random.seed(31371)


directory = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/nuc_gold/'
files  = os.listdir(directory)



translation_dict = dict()
codingsfilename = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/mark.txt'

with open(codingsfilename, 'r') as fin:
	for line in fin:
		copy = line[:-1]
		line_list = copy.split()
		translation_dict[line_list[1]] = list(map(int, line_list[2].split(',')))





import random
random.seed(7812)
randomwindow = [[0]*121, [0]*121]

def count_at(seq, pos, vec):
	for i in range(0, 121):
		e = pos - 60 + i
		if len(seq[e:e+3]) < 3:
			break
		else:
			vec[1][i] += 1
			if seq[e] in ['T', 'A']:
				vec[0][i] += 1
	return vec


matrix_st_no = list()
matrix_st_sh = list()
matrix_sh_sh = list()
matrix_start_sh = list()
matrix_start_no = list()


for filename in files:
	path = directory + filename
	bp_dict = reader.readfasta_todictionary(path)
	bp_dict = reader.filter_paralogues(bp_dict)

	for idr in bp_dict:
		tr = translation_dict[idr]
		seq = reader.clean(bp_dict[idr])

		if seq[-8:] == 'AAAAAAAA':
			while seq[-1] == 'A':
				seq = seq[:-1]

		if len(tr) > 2:
			matrix_st_sh.append(count_at(seq, tr[-1], [[0]*121, [0]*121]))
			matrix_start_sh.append( count_at(seq, tr[0], [[0]*121, [0]*121]))
			
			for shift in range(1, len(tr) - 1, 2):
				if tr[shift + 1] - tr[shift] == 1:
					matrix_sh_sh.append(count_at(seq, tr[shift], [[0]*121, [0]*121]))

		else:
			matrix_st_no.append(count_at(seq, tr[-1], [[0]*121, [0]*121]))
			matrix_start_no.append(count_at(seq, tr[0], [[0]*121, [0]*121]))

		# j += 1


def count_ratio(matrix, pos):
	up = 0 
	down = 0
	for elem in matrix:
		up += elem[0][pos]
		down += elem[1][pos]
	return up/float(down)

def count_ratios_prom(matrix, pos):
	bootstrap_list = list()
	for i in range(1000):
		
		newmatrix = list()
		for j in range(1000):
			newmatrix.append(random.choice(matrix))
		bootstrap_list.append(count_ratio(newmatrix, pos))
	bootstrap_list.sort()
	return bootstrap_list[49], count_ratio(matrix, pos), bootstrap_list[949]










def print_matrix(matrix, filename):
	with open(filename, 'w') as fout:
		for pos in range(121):
			low, mean, high = count_ratios_prom(matrix, pos)
			fout.write(str(pos - 60) + '\t' + str(low) + '\t' + str(mean)+ '\t' + str(high) + '\n' )


print_matrix(matrix_st_no, 'matrix_st_no_90.txt')
print_matrix(matrix_st_sh, 'matrix_st_sh_90.txt')
print_matrix(matrix_sh_sh, 'matrix_sh_sh_90.txt')
print_matrix(matrix_start_sh, 'matrix_start_sh_90.txt')
print_matrix(matrix_start_no, 'matrix_start_no_90.txt')


matrix_st_no = list()
matrix_st_sh = list()
matrix_sh_sh = list()
matrix_start_sh = list()
matrix_start_no = list()


# print_vec_sgl(st_no, 'at_st_no_sgl.txt')
# print_vec_sgl(st_sh, 'at_st_sh_sgl.txt')
# print_vec_sgl(sh_sh, 'at_sh_sh_sgl.txt')
# print_vec_sgl(start_no, 'at_start_no_sgl.txt')
# print_vec_sgl(start_sh, 'at_start_sh_sgl.txt')

# print_vec(st_no, 'at_st_no.txt')
# print_vec(st_sh, 'at_st_sh.txt')
# print_vec(sh_sh, 'at_sh_sh.txt')

