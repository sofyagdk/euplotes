import os 
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader

indirectory = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/nuc_gold/'
outdirectory = '/mnt/gamma/user/sofya/scripts/final/0pipeline/8_expression/animal_fastas/'

files = os.listdir(indirectory)
animals = ['petz', 'raik', 'octo', 'harp', 'eury', 'rari', 'foca', 'minu', 'cras']



analysed= 0 
for an in animals:
	outfile = open(outdirectory + an + '.fasta', 'w')
	for filename in files:
		path = indirectory + filename
		
		with open(path, 'r') as fin:
			E = fin.readlines()
			if len(E) == 0:
				continue
		
		bp_dict = reader.readfasta_todictionary(path)
		bp_dict = reader.filter_paralogues(bp_dict)
		for elem in bp_dict:
			bp_dict[elem] = reader.clean(bp_dict[elem])

		for idr in bp_dict:
			if idr[:4] == an:
				analysed += 1
				outfile.write('>' + idr + '\n')
				outfile.write(bp_dict[idr] + '\n')
	outfile.close()

print(analysed)