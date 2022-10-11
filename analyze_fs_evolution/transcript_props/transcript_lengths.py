import os
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader


directory = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/nuc_gold/'
files  = os.listdir(directory)
animals = [ 'raik', 'octo', 'harp', 'eury','rari', 'foca', 'minu', 'cras', 'petz']


LenDict = dict()
for an in animals:
	LenDict[an] = list()


for filename in files:
	path = directory + filename
	bp_dict = reader.readfasta_todictionary(path)

	for idr in bp_dict:
		seq = reader.clean(bp_dict[idr])
		LenDict[idr[:4]].append(len(seq))


for an in animals:
	print(an + ' ' + ','.join(list(map(str, LenDict[an]))))