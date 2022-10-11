import os
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader
import random 
from optparse import OptionParser
import numpy as np



codingsfilename = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/mark.txt'
translation_dict = dict()
with open(codingsfilename, 'r') as fin:
	for line in fin:
		copy = line[:-1]
		line_list = copy.split()
		translation_dict[line_list[1]] = list(map(int, line_list[2].split(',')))



directory = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/nuc_gold/'
files  = os.listdir(directory)

###################################################
## Functions
###################################################

def simulate_letter(Nprobs, pos):
	dot = random.random()
	passed = 0 
	for l in Nprobs:
		passed += Nprobs[l][pos]
		if passed > dot:
			break
	return l


def simulate_seq(Nprobs, length = 100):
	simulated_seq = str()
	
	for pos in range(length):
		simulated_seq +=  simulate_letter(Nprobs, pos)
	return simulated_seq


def count_stops(seq, size=None):
	if size == None:
		size = len(seq) - 2
	stops = [0]*size


	for pos in range(min(len(seq) - 2, size)):
		if seq[pos:pos+3] in ['TAA', 'TAG']:
			stops[pos] = (1)

	return np.array(stops)


###################################################
## Get true probabilities
###################################################


utr_list = list()

for filename in files:
	path = directory + filename
	bp_dict = reader.readfasta_todictionary(path)
	bp_dict = reader.filter_paralogues(bp_dict)

	for idr in bp_dict:

		seq = reader.clean(bp_dict[idr])
		coding = translation_dict[idr]
		utr = seq[coding[-1]+3:]

		if 'T' not in utr:
			if 'G' not in utr:
				if 'C' not in utr:
					continue


		if utr[-4:] == 'AAAA':
			PA = len(utr) - 4 
			while utr[PA - 1] == 'A':
				PA -= 1
			utr = utr[:PA]

		utr_list.append(utr)

Nprobs = { 
	'A': [0]*100,
	'T': [0]*100, 
	'G': [0]*100, 
	'C': [0]*100, 
}
denom = [0]*100



for elem in utr_list:
	for i in range(min(100, len(elem))):
		denom[i] += 1
		Nprobs[elem[i]][i] += 1

for N in Nprobs:
	for i in range(100):
		Nprobs[N][i] = Nprobs[N][i]/float(denom[i])



stop_probs = np.array([0]*98)
stop_denom = np.array([0]*98)

for elem in utr_list:

	stop_probs += count_stops(elem, size=98)
	denom = [1] * min(len(elem), 98)
	while len(denom) != 98:
		denom.append(0)
	stop_denom += np.array(denom)


# print(' '.join(list(map(str, stop_probs))))

# print(' '.join(list(map(str, np.true_divide(stop_probs,stop_denom)))))


# exit()

###################################################
## Bootstrap
###################################################


# def get_boostrap_set(some_list):
# 	boot_list = list()
# 	for i in range(len(some_list)):
# 		boot_list.append(random.choice(some_list))
# 	return boot_list


# for bootsrtap_rounds in range(1000):
# 	boostrap_set = get_boostrap_set(utr_list)

# 	for elem in boostrap_set:

# 		stop_probs += count_stops(elem, size=98)
# 		denom = [1] * min(len(elem), 98)
# 		while len(denom) != 98:
# 			denom.append(0)
# 		stop_denom += np.array(denom)


# 	print(' '.join(list(map(str, np.true_divide(stop_probs,stop_denom)))))

# exit()

###################################################
## Simulations
###################################################

Simulated_stops = [0]*100


for stop_probability_rounds in range(1000):
	stop_probs = np.array([0]*98)

	for sequnce_simulations_number in range(4000):
		simulated_seq = simulate_seq(Nprobs) 
		stop_probs += count_stops(simulated_seq)

	

	# print(' '.join(list(map(str, stop_probs))))
	stop_probs = stop_probs/float(4000)
	print(' '.join(list(map(str, stop_probs))))








# fout.close()