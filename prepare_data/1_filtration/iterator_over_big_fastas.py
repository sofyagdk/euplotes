
fasta  = 'concatenated_transcriptomes.fasta'




def fasta_iterator(fasta):
	with open(fasta, 'r') as fin:
		idr = str()
		seq = str()

		for line in fin:
			
			if '>' == line[0]:
				if idr:
					yield [idr, seq]
				seq = ''

				if '> ' == line[:2]:
					idr = line.strip().split()[1]
				else:
					idr = line.strip().split()[0][1:]
			else:
				seq += line.strip()


def ATComposition(seq):
	at = seq.count('A') + seq.count('a') + seq.count('T') + seq.count('t')
	return at/float(len(seq))

print('idr seq_len at_comp sp')
for elem in fasta_iterator(fasta):
	idr, seq = elem
	print(' '.join(map(str, [idr, len(seq), ATComposition(seq), idr[:4]])))



