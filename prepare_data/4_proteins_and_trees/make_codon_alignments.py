import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader
import os
import subprocess
import glob

###########################################
### Create dictionaries
###########################################

codings_true = dict()

with open('/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/mark.txt', 'r') as fin:
	for line in fin:
		cod = line.split()[2]
		cod = list(map(int, cod.split(',')))
		codings_true[line.split()[1]] = cod



###########################################
### Functions
###########################################


def write_fasta(seqdict, filename):
	with open(filename, 'w') as fout:
		for idr in seqdict:
			fout.write('>' + idr + '\n')
			fout.write(seqdict[idr] + '\n')


def translatorx(filename):

	line = ['perl', '/mnt/gamma/user/sofya/tools/translatorx.pl', '-i', filename,\
		 '-o', filename[:-6], '-p', 'M', '-t', 'F', '-w', '1', '-c', '10']	
	process = subprocess.Popen(line)
	process.communicate()


	translator_aligned = reader.readfasta_todictionary(filename[:-6] + ".nt_ali.fasta")
	

	# get a recursive list of file paths that matches pattern including sub directories
	fileList = glob.glob(filename[:-6] + "*")
	# Iterate over the list of filepaths & remove each file.
	for filePath in fileList:
	    try:
	        os.remove(filePath)
	    except OSError:
	        print("Error while deleting file")

	write_fasta(translator_aligned, filename)




###########################################
### MAIN
###########################################

already = os.listdir('./codon_alignments/')


directory = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/nuc_gold/'
files = os.listdir(directory)

for filename in files:
	if filename in already:
		continue
	path = directory + filename
	bp_dict = reader.readfasta_todictionary(path)
	codon_dict = dict()


	for idr in bp_dict:
		cod = codings_true[idr]
		seq = reader.clean(bp_dict[idr])
		print(cod)
		print(seq)

		codonedseq = str()
		for i in range(0, len(cod), 2):
			s = cod[i]
			e = cod[i+1]
			codonedseq  += seq[s:e]
		codon_dict[idr] = codonedseq



	write_fasta(codon_dict, './codon_alignments/' + filename)
	translatorx('./codon_alignments/' + filename)







