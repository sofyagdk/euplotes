

import os
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader


directory = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/nuc_gold/'
files  = os.listdir(directory)




# def fish_translated(idr):
# 	with open('/mnt/gamma/user/sofya/scripts/final/pictures/stops_fig2/translated_segs_8.txt') as fin:
# 		for line in fin:
# 			line_list = line.split()
# 			if line_list[0] == idr:

# 				cod = list(map(int, cod.split(',')))
# 				return cod
codingsfilename = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/mark.txt'
translation_dict = dict()
with open(codingsfilename, 'r') as fin:
	for line in fin:
		copy = line[:-1]
		line_list = copy.split()
		translation_dict[line_list[1]] = list(map(int, line_list[2].split(',')))


# hom_cras_c3120_g1_i1.fasta minu_c12482_g1_i1 20,827 20,827 40,848 40,848
no = list()
sh = list()



import random




for filename in files:
	path = directory + filename
	bp_dict = reader.readfasta_todictionary(path)
	bp_dict = reader.filter_paralogues(bp_dict)

	for idr in bp_dict:
		tr = translation_dict[idr]
		seq = reader.clean(bp_dict[idr])

		if seq[-7:] != 'AAAAAAA':
			continue
		PA = len(seq) - 7 
		while seq[PA - 1] == 'A':
			PA -= 1

		# print(seq)
		# print(PA)
		# pause = int(input())


		if len(tr) > 2:
			sh.append(PA - tr[-1])
			

		else:
			no.append(PA - tr[-1])




def print_vec(vec, filename):
	with open(filename, 'w') as fout:
		for elem in vec:
			fout.write(str(elem) + ' ' )

def print_vec_freq(vec, filename):
	with open(filename, 'w') as fout:
		for i in range(71):
			am = 0
			for elem in vec:
				if elem == i:
					am += 1
			fout.write(str(i) + ' ' + str( am/float(sum(vec) )) + '\n')


def print_vec_sgl(vec, filename):

	with open(filename, 'w') as fout:
		newvec = list()
		for i in range(71):
			am = 0
			for elem in vec:
				if elem == i:
					am += 1
			newvec.append(am/float(sum(vec)))

		fout.write(str(0) + ' ' + str( newvec[0]) + '\n')
		for i in range(1,70):
			fout.write(str(i) + ' ' + str( (newvec[i-1] + newvec[i]+ newvec[i+1])/3 ) + '\n')

		fout.write(str(70) + ' ' + str( newvec[70]) + '\n')

print_vec(no, 'PA_no_nums.txt')
print_vec(sh, 'PA_sh_nums.txt')

print_vec_freq(no, 'PA_no.txt')
print_vec_freq(sh, 'PA_sh.txt')
print_vec_sgl(no, 'PA_no_sgl.txt')
print_vec_sgl(sh, 'PA_sh_sgl.txt')
