#htis is the final variant

import sys
sys.path.append("../../../lib/")
import reader
import os
import numpy
import random
random.seed(2307)
from scipy.stats import chi2


def find_id(seq1, seq2):
	same = 0
	length = 0
	for i in range(len(seq1)):
		if seq1[i] != '-' and seq2[i] != '-':
			length += 1
			if seq1[i] == seq2[i]:
				same += 1
	return same, length


def big_sigma(data_list):
	n = len(data_list)
	variance = numpy.var(numpy.array(data_list))
	low, high = chi2.interval(0.95, n-1)
	big_sigma = float(variance *(n - 1)/low)
	return float(big_sigma**(0.5))




def cleanfromgaps(align_list):
	i = 0 
	while i < len(align_list[0]):
		if align_list[0][i] == '-':
			align_list[0] = align_list[0][:i] + align_list[0][i + 1:]
			align_list[1] = align_list[1][:i] + align_list[1][i + 1:]
		elif align_list[1][i] == '-':
			align_list[0] = align_list[0][:i] + align_list[0][i + 1:]
			align_list[1] = align_list[1][:i] + align_list[1][i + 1:]
		else:
			i += 1
	return align_list


def cut_to_ids(align_list, piece_len):
	data = list() 
	i = 0 
	for i in range(0, len(align_list[0]) - len(align_list[0])%piece_len, piece_len):
		a, b = find_id(align_list[0][i:i+piece_len], align_list[1][i:i+piece_len])
		data.append(a)
	return data



directory = '/mnt/mapr/user/mmoldovan/phosphosites/data/human_oggs_no_parals_mammalia_als/' ## must contain some protein alignments
files = os.listdir(directory)

starts = open('starts_vbigsigma.txt', 'w')


ends = open('ends_vbigsigma.txt', 'w')
outline = str()


# for piece_len in range(3, 61):
# 	outline += str(piece_len)+'m\t' + str(piece_len) + 'av\t' + str(piece_len) + 'sg\t'
# ends.write(outline[:-1] + '\n')
# starts.write(outline[:-1] + '\n')

ends.write('len m av sg\n')
starts.write('len m av sg\n')



for filename in files:
	path = directory + filename
	al_list = reader.readfasta_tolist(path)

	seq1, seq2 = random.sample(al_list, 2)
	seq1, seq2 = cleanfromgaps([seq1, seq2])

	
	for piece_len in range(3, min(61, len(seq1)//5)):
		my_meaning, b = find_id(seq1[-piece_len:], seq2[-piece_len:])
		endsdata = cut_to_ids([seq1[:-piece_len], seq2[:-piece_len]], piece_len)
		sigma = big_sigma(endsdata)
		av = numpy.mean(numpy.array(endsdata))

		# variance = numpy.var(numpy.array(endsdata))
		# sigma = variance ** (0.5)

		ends.write(' '.join(list(map(str, [piece_len, my_meaning, av, sigma])))  + '\n')


	outline = str()
	for piece_len in range(3, min(61, len(seq1)//5)):
		my_meaning, b = find_id(seq1[:piece_len], seq2[:piece_len])
		startsdata = cut_to_ids([seq1[piece_len: ], seq2[piece_len:]], piece_len)
		sigma = big_sigma(startsdata)
		av = numpy.mean(numpy.array(startsdata))

		# variance = numpy.var(numpy.array(startsdata))
		# sigma = variance ** (0.5)

		starts.write(' '.join(list(map(str, [piece_len, my_meaning, av, sigma])))  + '\n')

	






starts.close()
ends.close()


print('la finite')

