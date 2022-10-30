outfile = open('./ribobits.fasta', 'w')
num= 1
with open('/../../../data/E.crassus_Ribo.fq', 'r') as fin:
	for line in fin:
		if line[0] == '@':
			next = True
			continue
		elif (next):
			next = False
			outfile.write('>' + str(num) + '\n')
			outfile.write(line)
			num += 1
			continue
outfile.close()