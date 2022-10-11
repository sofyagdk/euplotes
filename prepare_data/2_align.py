import os
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader
import marking_module_v3 as mark
import bootstrap_tables_v4 as boot 



def makefasta_fromdict(bp_dict, fastaname):
	with open(fastaname, 'w') as fasta:
		for idr in bp_dict:
			fasta.write('>' + idr + '\n')
			fasta.write(bp_dict[idr] + '\n')

def align(fasta, alignmentfilename):
	cmd = 'muscle -in ' + fasta + ' -out ' + alignmentfilename #also folder name
	os.system(cmd)




# for filename in files:
	
# 	filename = indir + filename
# 	aligned = readfasta_todictionary(filename)

# 	# for idr in aligned:
# 	# 	seq = clean(aligned[idr])
# 	# 	coding, new_seq = Marking.getcoding(seq,  30, 30, 35)
# 	# 	aligned[idr] = new_seq
	
# 	with open('temporary_for_realign.fasta', 'w') as temp:
# 		for idr in aligned:
# 			temp.write('> ' + idr + '\n' )
# 			temp.write(str(aligned[idr]) + '\n' )
	
# 	align('temporary_for_realign.fasta', filename)



# ####the checker
# unrev_files = 0 
# for filename in files:
# 	unrev_local = 0 
# 	filename = indir + filename
# 	aligned = readfasta_todictionary(filename)
# 	for idr in aligned:
# 		seq = clean(aligned[idr])
# 		coding, new_seq = Marking.getcoding(seq,  30, 30, 35)
# 		if new_seq not in seq:
# 			unrev_local += 1
# 	if unrev_local != 0:
# 		unrev_files += 1
# 		# align(indir + filename, alignmentfile)

# print(unrev_files)

# def getreference_idr(idr): 
#     args = ["grep", str(idr), '/mnt/gamma/user/sofya/scripts/final/make_references/reference_pairs.txt']
#     line = subprocess.check_output(args)
#     line = line.decode("utf-8")

#     line = line.split('\t')
#     ref_idr = line[1]

#     return ref_idr



directory = './ortho_turned/'
files  = os.listdir(directory)

a = 0 
for filename in files:
	a += 1
	if a%100 == 0:
		print(a)

	path = directory + filename

	align(path, './nuc_alignments/' + filename)






