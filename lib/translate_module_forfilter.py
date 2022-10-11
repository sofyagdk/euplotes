#------------------------------------------------------------------------------------------------------
#TRANSLATE BLOCK
#------------------------------------------------------------------------------------------------------
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader
import math

translation_table = {
"TTT" : "F", "TCT" : "S", "TAT" : "Y", "TGT" : "C", 
"TTC" : "F", "TCC" : "S", "TAC" : "Y", "TGC" : "C", 
"TTA" : "L", "TCA" : "S", "TAA" : "*", "TGA" : "C", 
"TTG" : "L", "TCG" : "S", "TAG" : "*", "TGG" : "W", 
"CTT" : "L", "CCT" : "P", "CAT" : "H", "CGT" : "R", 
"CTC" : "L", "CCC" : "P", "CAC" : "H", "CGC" : "R", 
"CTA" : "L", "CCA" : "P", "CAA" : "Q", "CGA" : "R", 
"CTG" : "L", "CCG" : "P", "CAG" : "Q", "CGG" : "R", 
"ATT" : "I", "ACT" : "T", "AAT" : "N", "AGT" : "S", 
"ATC" : "I", "ACC" : "T", "AAC" : "N", "AGC" : "S", 
"ATA" : "I", "ACA" : "T", "AAA" : "K", "AGA" : "R", 
"ATG" : "M", "ACG" : "T", "AAG" : "K", "AGG" : "R", 
"GTT" : "V", "GCT" : "A", "GAT" : "D", "GGT" : "G", 
"GTC" : "V", "GCC" : "A", "GAC" : "D", "GGC" : "G", 
"GTA" : "V", "GCA" : "A", "GAA" : "E", "GGA" : "G", 
"GTG" : "V", "GCG" : "A", "GAG" : "E", "GGG" : "G" 
}

def translate(seq, coding):
	protein = str()
	for i in range(0, len(coding), 2):
		for codonstart in range(coding[i], coding[i+1], 3):
			protein += translation_table[seq[codonstart:codonstart+3]]
	return protein





def addtranslate(seq, tr_start, tr_end, protein):
	now = tr_start
	while now < 0:
		now += 3
	while now + 3 <= tr_end:
		protein += codons[seq[now:now+3]]
		now += 3
	if protein[0:4] == 'Stop':
		protein = protein[4:]
	return protein

def translatewithfs(seq, coding, idr):
	coding_s = list(coding)
	protein = str()
	idr_n = str(idr)
	if len(coding) == 0:
		return str(), idr_n
	# print('intranslate', coding)
	protein = addtranslate(seq, coding[0], coding[1], protein)

	while len(coding) > 2:
		idr_n = idr_n + '+shift+' + str(len(protein))
		coding = coding[2:]
		protein = addtranslate(seq, coding[0], coding[1], protein)

	# print(idr_n)
	return protein, idr_n

#------------------------------------------------------------------------------------------------------
#END OF TRANSLATE BLOCK
#------------------------------------------------------------------------------------------------------
