import os
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader


directory = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/nuc_gold/'
files  = os.listdir(directory)

codingsfilename = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/mark.txt'
truecodings_dict = dict()
with open(codingsfilename, 'r') as fin:
	for line in fin:
		copy = line[:-1]
		line_list = copy.split()
		truecodings_dict[line_list[1]] = list(map(int, line_list[2].split(',')))

c = 0 

animals = [ 'raik', 'octo', 'harp', 'eury','rari', 'foca', 'minu', 'cras']

# for an in animals:
# 	c = 0 
# 	apos = 0
# 	for elem in truecodings_dict:
# 		if elem[:4] != an:
# 			continue
# 		apos += 1
# 		if len(truecodings_dict[elem]) > 2:
# 			c += 1
# 	print(an)
# 	print(c)
# 	print(apos)
# 	# print(len(truecodings_dict))
# 	print(c/float(apos))
# # exit()




lys_contexts = list()


for filename in files:
	path = directory + filename
	bp_dict = reader.readfasta_todictionary(path)

	for idr in bp_dict:
		coding = truecodings_dict[idr]
		seq = reader.clean(bp_dict[idr])


		for i in range(0, len(coding), 2):
			piece = seq[coding[i]+3:coding[i+1]-3]
			for codon_start in range(12, len(piece) - 11, 3):
				codon = piece[codon_start:codon_start+3]
				if codon in ['AAA', 'AAG']:

					c = piece[codon_start - 9:codon_start + 12]
					lys_contexts.append(c)



with open('logo_lys.txt', 'w') as fout:
	for elem in lys_contexts:
		fout.write(elem + '\n')





context1 = list()
context2 = list()


for filename in files:
	path = directory + filename
	bp_dict = reader.readfasta_todictionary(path)

	for idr in bp_dict:
		coding = truecodings_dict[idr]
		seq = reader.clean(bp_dict[idr])

		if len(coding) < 3:
			continue

		for i in range(2, len(coding), 2):
			c = seq[coding[i-1] - 12:coding[i-1] + 10]

			if coding[i] - coding[i-1] == 2:
				context2.append(c)
				print(c, idr, i)
			else:
				context1.append(c)
				# print(c, idr, i)

# exit()

terminating = list()

for filename in files:
	path = directory + filename
	bp_dict = reader.readfasta_todictionary(path)

	for idr in bp_dict:
		coding = truecodings_dict[idr]
		seq = reader.clean(bp_dict[idr])


		if seq[-7:] != 'AAAAAAA':
			continue
		PA = len(seq) - 7 
		while seq[PA - 1] == 'A':
			PA -= 1

		seq = seq[:PA]

		c = seq[coding[-1] - 10:coding[-1] + 11]


		terminating.append(c)

with open('logo_1fs_tolys.txt', 'w') as fout:
	for elem in context1:
		if elem[9:12] not in ['AAA', 'AAG']:
			continue
		fout.write(elem[:12] + elem[13:] + '\n')


with open('logo_1fs.txt', 'w') as fout:
	for elem in context1:
		# if elem[9:12] not in ['AAA', 'AAG']:
		# 	continue
		fout.write(elem + '\n')

with open('logo_2fs.txt', 'w') as fout:
	for elem in context2:
		fout.write(elem + '\n')

with open('logo_term.txt', 'w') as fout:
	for elem in terminating:
		if len(elem) < 21:
			continue
		fout.write(elem + '\n')			