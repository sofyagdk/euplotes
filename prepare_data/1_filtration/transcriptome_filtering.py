import numpy as np
import pandas as pd


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

		yield [idr, seq]


def ATComposition(seq):
	at = seq.count('A') + seq.count('a') + seq.count('T') + seq.count('t')
	return at/float(len(seq))
	

fastas = ['E.euryhalinus_n.fasta', 
'E.focardii_n.fasta', 
'E.harpa_n.fasta', 
'E.minuta_n.fasta', 
'E.octocarinatus_n.fasta', 
'E.petzi_n.fasta', 
'E.raikovi_n.fasta', 
'E.rariseta_n.fasta', 
'T.thermophila_n.fasta', 
'E.crassus_combined.fasta']


# at_filter = {'petz': 0.4987458409862167,
#  'octo': 0.5967763899373861,
#  'raik': 0.5497183536720657,
#  'harp': 0.35742548155789683,
#  'minu': 0.5244974217432475,
#  'rari': 0.4747032781566095,
#  'cras': 0.5079885157345991,
#  'foca': 0.5816529135566422,
#  'eury': 0.4999457376578958,
#  'tetr': 0.6019639117959017}


# read_filter_table:
# that has -- identificator + okay or not 

# filter_table = pd.read_csv('filter_table.txt', header = 0, sep = ' ')

# def IsOkay(idr):
# 	temp = filter_table[filter_table['idr'] == idr]
# 	return list(temp['any'])[0]


# for sp_fasta in fastas:
# 	new_fasta = open('./filtered/f'+sp_fasta, 'w') 


# 	for elem in fasta_iterator('./primary_files/'+sp_fasta):
# 		idr, seq = elem

# 		if IsOkay(idr):
# 			new_fasta.write('>' + idr + '\n')
# 			new_fasta.write(seq + '\n')

# 	new_fasta.close()




print('idr seq_len at_comp sp')

for sp_fasta in fastas:
	for elem in fasta_iterator('./filtered/f'+ sp_fasta):
		idr, seq = elem
		sp = idr[:4]
		at_comp = ATComposition(seq)
		out = [idr, len(seq), at_comp, sp]
		print(' '.join(map(str, out)))




