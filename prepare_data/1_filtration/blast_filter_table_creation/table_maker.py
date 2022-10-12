#this is the file that allows to analyze chimeras! 

import os 

#----------------------------------------------------------------
# creating genuses list 
#----------------------------------------------------------------

def getgenus(filename):
	if '_' in filename:
		genus = filename.split('_')
		genus = genus[0]
	else:
		genus = filename[:13]
	return genus


directory = '/mnt/gamma/user/sofya/data/RNA_for_blast/tsvs_id70_e-25/'
dbfiles = os.listdir(directory)
genuses = list()
for filename in dbfiles:
	genus = getgenus(filename)
	if genus not in genuses:
		genuses.append(genus)


#----------------------------------------------------------------
# creating list of transcriptes
#----------------------------------------------------------------

idrs = dict()

with open('/mnt/gamma/user/sofya/scripts/comb/cilia_having_ortholouges.fasta', 'r') as fin:
	for line in fin:
		if line[0] == '>':
			idr = line[1:]
			idr = idr.split()
			idr = idr[0]
			idrs[idr] = list()

#----------------------------------------------------------------
# collecting the data from tsvs to dictionary 
#----------------------------------------------------------------


for adb in dbfiles:
	with open(directory + adb, 'r') as fin:
		genus = getgenus(adb)
		
		for line in fin:
			if line[0] == '#':
				continue
			llist = line.split()
			query = llist[0]
			qstart = llist[6]
			qend = llist[7]
			identity =  float(llist[2])
			idrs[query].append(tuple([genus, identity, qstart, qend]))

# #----------------------------------------------------------------
# # printing the dictionary. LENGTHS
# #----------------------------------------------------------------

# print(genuses)

# with open('cilia_vs_animals_70_e-25_v0.txt', 'w') as fout:
# 	fout.write('idr\t')
# 	outline = 'idr\t'
# 	for genus in genuses:
# 		outline += genus + '_s\t' + genus + '_e\t'
# 	outline += 'animal\n'
# 	fout.write(outline) 


# 	for idr in idrs:
# 		outline = idr + '\t'

# 		for genus in genuses:
# 			biggestid = 0 
# 			st = 0 
# 			en = 0

# 			for tup in idrs[idr]:

# 				if tup[0] == genus:
# 					if tup[1] > biggestid:
# 						biggestid = tup[1]
# 						st = tup[2]
# 						en = tup[3]
			
# 			outline += str(st) + '\t' + str(en) + '\t'
		
# 		outline += idr[:4]
# 		outline += '\n'
# 		fout.write(outline)

#----------------------------------------------------------------
# printing the dictionary. IDENTITY
#----------------------------------------------------------------

print(genuses)

with open('cilia_vs_animals_70_e-25_v1_identity.txt', 'w') as fout:
	fout.write('idr\t')
	outline = 'idr\t'
	for genus in genuses:
		outline += genus + '\t'
	outline += 'animal\n'
	fout.write(outline) 


	for idr in idrs:
		outline = idr + '\t'

		for genus in genuses:
			biggestid = 0 

			for tup in idrs[idr]:

				if tup[0] == genus:
					if tup[1] > biggestid:
						biggestid = tup[1]
			
			outline += str(biggestid) + '\t'
		
		outline += idr[:4]
		outline += '\n'
		fout.write(outline)

