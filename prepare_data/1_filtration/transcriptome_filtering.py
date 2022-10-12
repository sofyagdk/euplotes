####################################
# sofya 08.04.2022
# sofya.gaydukova@gmail.com

# transcriptome filter
# uses filter table and pre-calculated at composition thresholds
####################################

import pandas as pd
import sys
sys.path.append('../../lib/')
from reader import fasta_iterator


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


## values calculated from AT distributions of transcriptomes
at_filter = {'petz': 0.4987458409862167,
 'octo': 0.5967763899373861,
 'raik': 0.5497183536720657,
 'harp': 0.35742548155789683,
 'minu': 0.5244974217432475,
 'rari': 0.4747032781566095,
 'cras': 0.5079885157345991,
 'foca': 0.5816529135566422,
 'eury': 0.4999457376578958,
 'tetr': 0.6019639117959017}


filter_table = pd.read_csv('filter_table.txt', header=0, sep=' ')

def IsOkay(idr):
	temp = filter_table[filter_table['idr'] == idr]
	return list(temp['any'])[0]


for sp_fasta in fastas:
	new_fasta = open('./filtered/f'+sp_fasta, 'w')


	for elem in fasta_iterator('../data/primary_files/'+sp_fasta):
		idr, seq = elem

		if IsOkay(idr):
			new_fasta.write('>' + idr + '\n')
			new_fasta.write(seq + '\n')

	new_fasta.close()
