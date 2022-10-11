import os
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader
import random 
from optparse import OptionParser




#######################################################################################
########## Data tables
#######################################################################################
translation_table = {
"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
"TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
"TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
"TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V", 
"TCT" : "S", "CCT" : "P", "ACT" : "T", "GCT" : "A",
"TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
"TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
"TCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
"TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
"TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
"TAA" : "*", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
"TAG" : "*", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
"TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G", 
"TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
"TGA" : "C", #u-hoo important here
"CGA" : "R", "AGA" : "R", "GGA" : "G",
"TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G"
}

codon_list = [
"TTT", "TCT", "TAT", "TGT", 
"TTC", "TCC", "TAC", "TGC", 
"TTA", "TCA", "TAA", "TGA", 
"TTG", "TCG", "TAG", "TGG", 
"CTT", "CCT", "CAT", "CGT", 
"CTC", "CCC", "CAC", "CGC", 
"CTA", "CCA", "CAA", "CGA", 
"CTG", "CCG", "CAG", "CGG", 
"ATT", "ACT", "AAT", "AGT", 
"ATC", "ACC", "AAC", "AGC", 
"ATA", "ACA", "AAA", "AGA", 
"ATG", "ACG", "AAG", "AGG", 
"GTT", "GCT", "GAT", "GGT", 
"GTC", "GCC", "GAC", "GGC", 
"GTA", "GCA", "GAA", "GGA", 
"GTG", "GCG", "GAG", "GGG" 
]

aminoacid_code = {'Y': ['TAT', 'TAC'], 'C': ['TGC', 'TGT', 'TGA'], 
'F': ['TTT', 'TTC'], 'I': ['ATC', 'ATT', 'ATA'], 
'L': ['CTC', 'CTT', 'CTG', 'CTA', 'TTG', 'TTA'], 
'K': ['AAG', 'AAA'], 'H': ['CAC', 'CAT'], 'G': ['GGT', 'GGC', 'GGG', 'GGA'], 
'D': ['GAC', 'GAT'], 'V': ['GTT', 'GTG', 'GTA', 'GTC'], 'X': ['TAG', 'TAA'], 
'P': ['CCC', 'CCG', 'CCA', 'CCT'], 'A': ['GCC', 'GCG', 'GCA', 'GCT'], 
'W': ['TGG'], 'M': ['ATG'], 'Q': ['CAA', 'CAG'],
'R': ['AGA', 'CGG', 'CGC', 'CGT', 'AGG', 'CGA'], 
'S': ['TCC', 'TCT', 'TCG', 'AGT', 'AGC', 'TCA'], 'N': ['AAC', 'AAT'], 
'E': ['GAA', 'GAG'], 'T': ['ACT', 'ACC', 'ACG', 'ACA']}

aminoacids_list = ['Y', 'C', 'F', 'I', 'L', 
'K', 'H', 'G', 'D', 'V', 
'P', 'A', 'W', 'M', 'Q', 
'R', 'S', 'N', 'E', 'T']


bp = ['A', 'T', 'G', 'C']

species = ['cras', 'minu', 'foca', 'rari', 
		'eury', 'harp', 'octo', 'raik'] # NO PETZ!




#######################################################################################
########## Functions
#######################################################################################

def get_refseqset(expr_dict, Sdict, procent, sp):
	expression_pairs = list()
	for idr in expr_dict:
		if idr[:4] != sp:
			continue
		expression_pairs.append((idr, expr_dict[idr]))

	expression_pairs.sort(key = lambda x: x[1], reverse = True)
	refseqset = list()
	refseqidrs = list()
	for i in range(0, int(len(expression_pairs)*procent)):
		refseqset.append(Sdict[expression_pairs[i][0]])
		refseqidrs.append(expression_pairs[i][0])

	return refseqidrs, refseqset




def get_codon_relative_adaptiveness(refseqset):
	Fdict = dict()
	for codon in codon_list:
		Fdict[codon] = 0
	codon_num = 0

	for seq in refseqset:
		codon_num += len(seq)//3

		for i in range(0, len(seq), 3):
			codon = seq[i:i+3]
			Fdict[codon] += 1

	for codon in codon_list:
		Fdict[codon] /= float(codon_num)


	codon_relative_adaptiveness = dict()

	for aa in aminoacids_list:
		aa_codon_freqs = list()

		for aa_codon in aminoacid_code[aa]:
			aa_codon_freqs.append(Fdict[aa_codon])
		most_frequent = max(aa_codon_freqs)

		for aa_codon in aminoacid_code[aa]:
			if most_frequent == 0:
				continue
			codon_relative_adaptiveness[aa_codon] = Fdict[aa_codon]/most_frequent

	return codon_relative_adaptiveness


def calculate_CAI(gene_seq, codon_relative_adaptiveness):
	L = len(gene_seq)//3
	CAI = 1
	for i in range(0, len(gene_seq), 3):
		codon = gene_seq[i:i+3]
		CAI *= float(codon_relative_adaptiveness[codon])
		
	CAI = CAI ** (1/float(L))

	return CAI

# 1.6711302142e-51
# 0.751888917617





#######################################################################################
########## get all dictionaries
#######################################################################################

codingsfilename = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/mark.txt'

def get_norm_codings_dict(codingsfilename):
    nco_dict = dict()
    with open(codingsfilename, 'r') as fin:
        for line in fin:
            copy = line.split()
            #hom_cras_c9739_g1_i1.fasta minu_c24042_g1_i1 87,402,404,1307 272,1307 87,402,404,1307 272,1307
            coding = copy[2].split(',')
            coding =  list(map(int, coding))
            nco_dict[copy[1]] = coding
    return nco_dict


nco_dict = get_norm_codings_dict(codingsfilename)

def get_codoned_seq(seq, coding):
	cseq = str()
	for i in range(0,len(coding),2):
		cseq = seq[coding[i]:coding[i+1]]
	return cseq


expr_dict = dict() #idr - tpm
with open('/mnt/gamma/user/sofya/scripts/final/0pipeline/8_expression/abundance_all.tsv', 'r') as fin:
	for line in fin:
		if line[:6] == 'target':
			continue
		copy = line.split()
		expr_dict[copy[0]] = float(copy[4])


exp_table = dict() #idr animal length eff_length est_counts tpm coding main coding_ali main_ali fs pnps nd_pnps N_pnps sd_pnps S_pnps w dN N dS S
with open('/mnt/gamma/user/sofya/scripts/final/0pipeline/8_expression/expt/expression_table.txt', 'r') as fin:
	next(fin)
	for line in fin:
		if line[:3] == 'idr':
			continue
		copy = line.split()
		exp_table[copy[0]] = line.strip()




directory = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/nuc_gold/'
files  = os.listdir(directory)

Sdict = dict() #identificator - coding sequence (by codons)

for filename in files:
# for filename in files:

	path = directory + filename
	bp_dict = reader.readfasta_todictionary(path)
	bp_dict = reader.filter_paralogues(bp_dict)
	for idr in bp_dict:
		seq = reader.clean(bp_dict[idr])
		coding = nco_dict[idr]
		cseq = get_codoned_seq(seq, coding)
		Sdict[idr] = cseq 










#######################################################################################
########## write results
#######################################################################################

fout = open('expression_table_cai.txt', 'w')
fout.write('idr animal length eff_length est_counts tpm coding main coding_ali main_ali fs pnps nd_pnps N_pnps sd_pnps S_pnps w dN N dS S 5refsp 5sp 10refsp 10sp 20refsp 20sp 5 10 20\n')


addwrite = dict()
for idr in Sdict:
	addwrite[idr] = list()

for sp in species:
	for procent in [0.05, 0.1, 0.2]:
		refseqidrs, refseqset = get_refseqset(expr_dict, Sdict, procent, sp)
		codon_relative_adaptiveness = get_codon_relative_adaptiveness(refseqset)
		for idr in Sdict:
			if idr[:4] == sp:
				cai = calculate_CAI(Sdict[idr], codon_relative_adaptiveness)
				addwrite[idr].append(int(idr in refseqidrs))
				addwrite[idr].append(cai)



for procent in [0.05, 0.1, 0.2]:
	sumrefseqset = list()
	for sp in species:
		refseqidrs, refseqset = get_refseqset(expr_dict, Sdict, procent, sp)
		sumrefseqset.extend(refseqset)
	codon_relative_adaptiveness = get_codon_relative_adaptiveness(sumrefseqset)
	for idr in Sdict:
		cai = calculate_CAI(Sdict[idr], codon_relative_adaptiveness)
		addwrite[idr].append(cai)




for idr in Sdict:
	fout.write(exp_table[idr])
	fout.write(' ')
	if idr[:4] == 'petz':
		fout.write('      ')	
		fout.write(' '.join(list(map(str, addwrite[idr]))))
	else:
		fout.write(' '.join(list(map(str, addwrite[idr])))) 
	fout.write('\n')

