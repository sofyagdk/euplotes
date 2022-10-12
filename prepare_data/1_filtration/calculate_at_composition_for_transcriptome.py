####################################
# sofya 08.04.2022
# sofya.gaydukova@gmail.com

# get information about transcriptomes
# print into file at_comp_after_filter.txt or at_comp_before_filter.txt
####################################

import sys
sys.path.append('../../lib/')
from reader import fasta_iterator


def ATComposition(seq):
	at = seq.count('A') + seq.count('a') + seq.count('T') + seq.count('t')
	return at/float(len(seq))

fasta  = '../data/concatenated_transcriptomes.fasta'

print('idr seq_len at_comp sp')
for elem in fasta_iterator(fasta):
	idr, seq = elem
	print(' '.join(map(str, [idr, len(seq), ATComposition(seq), idr[:4]])))
