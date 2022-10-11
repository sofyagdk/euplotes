####################################
# sofya 08.04.2022
# sofya.gaydukova@gmail.com

# reads files from ./harpa_para/
# containing paralogs
# turns and alignes files 
# stores into folder ./harpa_para_1turned/
####################################


import os
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader
import marking_module_v3 as mark
import bootstrap_tables_v4 as boot 


####################################
# FUNCTIONS
####################################



def makefasta_fromdict(bp_dict, fastaname):
	with open(fastaname, 'w') as fasta:
		for idr in bp_dict:
			fasta.write('>' + idr + '\n')
			fasta.write(bp_dict[idr] + '\n')

def align(fasta, alignmentfilename):
	cmd = 'muscle -in ' + fasta + ' -out ' + alignmentfilename #also folder name
	os.system(cmd)


comp_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
def ReverseComplimentary(seq):
	cseq = ''
	for bp in seq:
		cseq += comp_dict[bp]
	return cseq[::-1]



def TurnOnePara(para_dict):
	idrs = para_dict.keys()
	if para_dict[idrs[0]][-6:] == 'AAAAAA':
		para_dict[idrs[1]]  = ReverseComplimentary(para_dict[idrs[1]])
	elif para_dict[idrs[0]][:6] == 'TTTTTT':
		para_dict[idrs[0]]  = ReverseComplimentary(para_dict[idrs[0]])
	elif para_dict[idrs[1]][:6] == 'TTTTTT':
		para_dict[idrs[1]]  = ReverseComplimentary(para_dict[idrs[1]])
	elif para_dict[idrs[1]][-6:] == 'AAAAAA':
		para_dict[idrs[0]]  = ReverseComplimentary(para_dict[idrs[0]])
	else:
		para_dict[idrs[0]]  = ReverseComplimentary(para_dict[idrs[0]])
	return para_dict



####################################
# FUNCTIONS
####################################



directory = './harpa_para/'
filenames  = os.listdir(directory)

outdirectory = './harpa_para_1turned/'
a = 0 
for filename in filenames:

	path = directory + filename
	bp_dict = reader.readfasta_todictionary(path)
	bp_dict  = TurnOnePara(bp_dict )

	makefasta_fromdict(bp_dict, outdirectory + filename)


	align(outdirectory + filename, outdirectory + filename)



