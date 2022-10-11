def readfasta_todictionary(filename): # changed 19-11-27 so could read water fasta output
	seqdict = dict()
	fastafile = open(filename, 'r')
	idr = 'identificator'
	seq = 'sequence'
	for line in fastafile:
		if line[0] == '#':
			continue
		if line[0] == '>':
			seqdict[idr] = seq
			if line[:2] == '> ':
				idr = line.split()[1]
			else:
				idr = line[1: -1].split()[0]
			# try: 
			# 	idr = line.split()[1]
			# except IndexError:
			# 	idr = line[1: -1]
			seq = str()
		else:
			seq += line[:-1]
	del seqdict['identificator']
	seqdict[idr] = seq
	fastafile.close()
	return seqdict


def writefasta_fromdictionary(fasta_dict, filename):
	with open(filename, 'w') as fout:
		for idr in fasta_dict:
			fout.write('> '+ idr + '\n')
			fout.write(fasta_dict[idr] + '\n')


def readfasta_tolist(filename): # changed 19-11-27 so could read water fasta output
	seqlist = list()
	fastafile = open(filename, 'r')
	currseq = str()

	with open(filename, 'r') as fastafile:
		for line in fastafile:
			if line[0] == '>':
				if currseq == '':
					continue
				else:
					seqlist.append(currseq)
					currseq = str()
			else:
				currseq += line[:-1]
	seqlist.append(currseq)

	return seqlist


def filter_paralogues(bp_dict):
	animals = list()
	para = list()
	new_dict = dict()

	for key in bp_dict.keys():
		animal  = key[:4]
		if animal in animals:
			para.append(animal)
		animals.append(key[:4])

	for key in bp_dict:
		animal  = key[:4]
		if animal in para:
			continue
		else:
			new_dict[key] = bp_dict[key]

	return new_dict

def choose_reference(bp_dict):
	ref = ''
	ref_len = 0
	for idr in bp_dict:
		if len(clean(bp_dict[idr])) > ref_len:
			ref = idr
			ref_len = len(clean(bp_dict[idr]))
	return ref

def clean(seq):
	new_seq = str()
	for elem in seq:
		if elem != '-':
			new_seq += elem
	return new_seq

def complementary(seq):
	newseq = str()
	seq = seq[::-1]
	compl = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', 'M':'K','R':'Y', 'W':'W', 'S':'S', 'Y':'R','K':'M', 'V':'B', 'H':'D', 'D':'H', 'B':'V'}
	for i in seq:
		newseq += compl[i]
	return newseq
# def getcoding(seqid,  start, back, forw):
# 	filename = ''.join(['/mnt/gamma/user/sofya/scripts/final/seq_tables/filter/', seqid, '.txt'])
# 	re = '.'.join(list(map(str, [start, back, forw])))	
# 	process = subprocess.Popen(["grep", re, filename], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
# 	stdout, stderr = process.communicate()
# 	stdout = stdout.decode('utf-8')
# 	# print(stdout)
# 	coding_line = stdout.split('[')[2]
# 	coding_line = coding_line.split(']')[0]
# 	coding_events = list(map(int, coding_line.split(',')))
# 	return coding_events