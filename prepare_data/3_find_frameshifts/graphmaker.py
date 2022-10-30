
### different from the previous variant only with the dictionaries!

# this needed for the algorithmic part
import random
import os
import glob
import subprocess
import numpy as np

# this needed for the main part
import os
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader
import random 
from optparse import OptionParser
import ds



def clean(seq):
	new_seq = str()
	for elem in seq:
		if elem != '-':
			new_seq += elem
	return new_seq
codons = {
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

def cut_to_ids(align_list, piece_len):
	data = list() 
	i = 0 
	for i in range(0, len(align_list[0]) - len(align_list[0])%piece_len, piece_len):
		a, b = find_id(align_list[0][i:i+piece_len], align_list[1][i:i+piece_len])
		if b == 0:
			continue
		data.append(float(a)/float(b))
	return data


after_sigmas = {4: 3.5574387539990515,
 5: 3.6213616483955673,
 6: 3.6839247473616443,
 7: 3.7671966511899,
 8: 3.860994266059596,
 9: 3.9313657054657334,
 10: 3.996616746642646,
 11: 4.02240567695236,
 12: 4.071695276167881,
 13: 4.086989989947383,
 14: 4.111182425510172,
 15: 4.143666419568383,
 16: 4.176671918620213,
 17: 4.227696667025651,
 18: 4.235057693533199,
 19: 4.270658503109524,
 20: 4.275957267967549,
 21: 4.298300097634783,
 22: 4.3010967838096406,
 23: 4.328965980064519,
 24: 4.36090515070121,
 25: 4.409000530122587,
 26: 4.425913603241926,
 27: 4.45741804056181,
 28: 4.47547432444811,
 29: 4.500000000000002,
 30: 4.48701419334597,
 31: 4.502418548794046,
 32: 4.485761787973615,
 33: 4.49551413538128,
 34: 4.486737785912649,
 35: 4.520672009161797,
 36: 4.537361244076282,
 37: 4.573004557626885,
 38: 4.58737116691905,
 39: 4.623108729728866,
 40: 4.636730652848285,
 41: 4.644687298967185,
 42: 4.658965880431205,
 43: 4.724934155485667,
 44: 4.776045345085099,
 45: 4.794948969396977,
 46: 4.765221455430356,
 47: 4.763613863400693,
 48: 4.757787476760424,
 49: 4.753155060528982,
 50: 4.757688890896863,
 51: 4.752790289256523,
 52: 4.7385764282668354,
 53: 4.76679980018736,
 54: 4.776747307193779,
 55: 4.787648305611946,
 56: 4.733472961694012,
 57: 4.726919545569095,
 58: 4.714381376607965,
 59: 4.736479824593554,
 60: 4.769353865296826}

before_sigmas = {4: 3.915879529338936,
 5: 4.045776276268301,
 6: 4.158067808186772,
 7: 4.263518536875723,
 8: 4.395319846136616,
 9: 4.495857845081842,
 10: 4.578546742640129,
 11: 4.666281292248686,
 12: 4.759369421308871,
 13: 4.871388069963927,
 14: 4.962157675795862,
 15: 5.044561379948813,
 16: 5.11395021593293,
 17: 5.193646027724873,
 18: 5.283888509926305,
 19: 5.401555967487762,
 20: 5.462897303465483,
 21: 5.532833804658966,
 22: 5.575833112172046,
 23: 5.653059481241112,
 24: 5.731187294079608,
 25: 5.805426858523965,
 26: 5.9005590795839415,
 27: 5.953871747482821,
 28: 5.996789951303679,
 29: 6.0352555521666496,
 30: 6.085258756100974,
 31: 6.12523582809225,
 32: 6.132997230458887,
 33: 6.102453241794649,
 34: 6.142637301245863,
 35: 6.178296943901003,
 36: 6.277039085837626,
 37: 6.311053089975982,
 38: 6.37252202938398,
 39: 6.38233419507839,
 40: 6.411512286408133,
 41: 6.392289559601434,
 42: 6.411127723593137,
 43: 6.376109756344804,
 44: 6.457831567891522,
 45: 6.545323887704407,
 46: 6.639799811857201,
 47: 6.595588765427529,
 48: 6.558424612318151,
 49: 6.530161792864341,
 50: 6.512096810160585,
 51: 6.499445081512664,
 52: 6.551602472553387,
 53: 6.668913813877647,
 54: 6.783695948229052,
 55: 6.8704784379359625,
 56: 6.875525636693364,
 57: 6.895065425253189,
 58: 6.85454914710484,
 59: 6.843490154530841,
 60: 6.772962119567572}

def create_filename(flag = None, extention = '.txt', num = 15):
		filename = str(flag)
		for i in range(num):
			filename += random.choice('0123456789') 
		filename += extention
		return filename


def translate_frame(seq):
	protein = str()
	for i in range(0, len(seq), 3):
		protein += codons[seq[i:i+3]]
	if "*" in protein:
		print("WARNING::graphlike::translate_frame:: the sequence contained a STOP codon.")
	return protein

def find_id(seq1, seq2):
	same = 0
	ungapped = 0
	for i in range(len(seq1)):
		if seq1[i] != '-' and seq2[i] != '-':
			ungapped += 1
			if seq1[i] == seq2[i]:
				same += 1
	return same, ungapped

# def align_local(fasta, idr_list, alignmentfilename):
# 	inline = str()
# 	for idr in idr_list:
# 		inline += fasta + ':' + idr + ' '

# 	cmd = 'water ' + inline + ' -outfile ' + alignmentfilename + ' -gapopen 10 -gapextend 0.5 -aformat3 fasta -auto'#also folder name

# 	# cmd = "muscle", "-in", filename, "-out", filename, "-quiet"
# 	os.system(cmd)



def get_p(nd1, N1, sd1, S1,nd2, N2, sd2, S2):
    nd1_ = np.random.poisson(lam = nd1, size = 1000)
    N1_ = np.random.poisson(lam = N1, size = 1000)
    sd1_ = np.random.poisson(lam = sd1, size = 1000)
    S1_ = np.random.poisson(lam = S1, size = 1000)
    nd2_ = np.random.poisson(lam = nd2, size = 1000)
    N2_ = np.random.poisson(lam = N2, size = 1000)
    sd2_ = np.random.poisson(lam = sd2, size = 1000)
    S2_ = np.random.poisson(lam = S2, size = 1000)
    
    sd1_ = np.where(sd1_==0, 0.001, sd1_) 
    sd2_ = np.where(sd2_==0, 0.001, sd2_) 
    N1_ = np.where(N1_==0, 0.001, N1_) 
    N2_ = np.where(N2_==0, 0.001, N2_) 

    pnps1 = nd1_/N1_/sd1_*S1_
    pnps2 = nd2_/N2_/sd2_*S2_
    procent_of_second_bigger = float(sum(pnps2>pnps1))/float(1000)
    return procent_of_second_bigger
    





class Frame(object):

	def __init__(self, start, end, parent, tomain = 0, is3 = False):
		self._start = start# coordinate of the first Nucleotide in a nuc alignment
		self._end = end # coordinate of the last + 1 Nucleotide in a nuc alignment
		self._parent  = parent# a nucleotide alignment sequence
		self._alignment_nuc = parent[self._start:self._end]
		self._arrows = list()
		self.out = dict()
		self._tomain = tomain
		self._is3 = is3
		self.out['_start'] = self._start
		self.out['_end'] =  self._end
		self.out['is3']= is3
		self.out['tomain'] = tomain



	def seq_nuc(self):
		return clean(self._alignment_nuc)

	def seq_prot(self):
		return  translate_frame(clean(self._alignment_nuc))

	def add_arrow(self, OtherFrame):
		self._arrows.append(OtherFrame)

	def get_length(self):
		l = len(clean(self._alignment_nuc))
		self.out['length'] = l
		return l

	def nuc_id(parent_pair):
		frame_pair = parent_pair[self._start:self._end]
		same, ungapped = find_id(self._alignment_nuc, frame_pair)
		if ungapped == 0:
			print("WARNING::Frame class::nuc_id::  frame has 0 ungapped positions in the nucleotide alignment")
			self.out['nuc_id'] = [0, 0]
			return 0
		if ungapped < 60:
			self.out['nuc_id']
			print("WARNING::Frame class::nuc_id::  frame has {0} ungapped nucleotide positions in the nucleotide alignment".format(str(ungapped)))
		self.out['nuc_id'] = [same, ungapped]
		return same/float(ungapped)



	def find_pair_protein_id(self, OtherFrame):
		filename = create_filename(flag = "frame_prot_al_", extention = '.fasta')
		
		with open(filename, 'w') as temp:
			temp.write("> native " + str(self._start) + 'to' + str(self._end) + '\n')
			temp.write(self.seq_prot() + '\n')
			temp.write("> pair " + str(OtherFrame._start) + 'to' + str(OtherFrame._end) + '\n')
			temp.write(OtherFrame.seq_prot() + '\n')
		
		# cmd = "muscle", "-in", filename, "-out", filename, "-quiet -stable"

		inline = str()
		inline += filename + ':native ' + filename + ':pair '
		print(inline)
		cmd = 'water ' + inline + ' -outfile ' + filename + ' -gapopen 10 -gapextend 0.5 -aformat3 fasta -sprotein'#also folder name
		print(cmd)


		os.system(cmd)
		water_aligned = reader.readfasta_tolist(filename)
		self.last_water_protein_alignment = water_aligned
		os.remove(filename)
		
		same, ungapped = find_id(water_aligned[0], water_aligned[1])
		if (('prot_ids' in self.out) == False):
			self.out['prot_ids'] = list()

		if ungapped == 0:
			print("WARNING::Frame class::find_pair_protein_id:: frame (:  has 0 ungapped positions in the protein alignment")
			# self.outall.append(0, 0)
			return 0, 0
		if ungapped < 20:
			print("WARNING::Frame class::find_pair_protein_id:: frame (:  has {0} ungapped aminoacid positions in the protein alignment".format(str(ungapped)))
		
		if (([same, ungapped] in self.out['prot_ids']) == False):
			self.out['prot_ids'].append([same, ungapped])
		return same, ungapped

	# def find_pair_protein_id_difference(self, OtherFrame):
	# 	my_meaning, b = boot.find_id(seq1[-piece_len:], seq2[-piece_len:])
	# 	if b == 0:
	# 		continue
	# 	my_meaning = float(my_meaning)/float(b)
	# 	endsdata = cut_to_ids([seq1[:-piece_len], seq2[:-piece_len]], piece_len)
	# 	sigma = boot.big_sigma(endsdata)
	# 	av = numpy.mean(numpy.array(endsdata))

	# 	variance = numpy.var(numpy.array(endsdata))
	# 	sigma = variance ** (0.5)


		





	def find_pair_pnps(self, OtherFrame):
		filename = create_filename(flag = "frame_pnps_al_", extention = '.fasta')

		with open(filename, 'w') as temp:
			temp.write("> native " + str(self._start) + 'to' + str(self._end) + '\n')
			temp.write(self.seq_nuc() + '\n')
			temp.write("> pair " + str(OtherFrame._start) + 'to' + str(OtherFrame._end) + '\n')
			temp.write(OtherFrame.seq_nuc() + '\n')

		
		line = ['perl', '/mnt/gamma/user/sofya/tools/translatorx.pl', '-i', filename,\
		 '-o', filename[:-6], '-p', 'M', '-t', 'F', '-w', '1', '-c', '10']	
		process = subprocess.Popen(line)
		process.communicate()


		translator_aligned = reader.readfasta_tolist(filename[:-6] + ".nt_ali.fasta")
		

		# get a recursive list of file paths that matches pattern including sub directories
		fileList = glob.glob(filename[:-6] + "*")
		# Iterate over the list of filepaths & remove each file.
		for filePath in fileList:
		    try:
		        os.remove(filePath)
		    except OSError:
		        print("Error while deleting file")
		# os.remove(filename[:-6] + "*")

		
		nd, N, sd, S = ds.getdnds(translator_aligned[0], translator_aligned[1])
		try: 
			pnps = float(nd)/ float(N) / float(sd) * float(S)
		except ZeroDivisionError:
			pnps = None

		if (('dnds' in self.out) == False):
			self.out['dnds'] = list()
		self.out['dnds'].append([pnps, nd, N, sd, S])
		
		return pnps, [nd, N, sd, S]   #maybe we will need to return the whole



	def find_paml_w(self, OtherFrame):

		filename = create_filename(flag = "paml_w_", extention = '.fasta')

		with open(filename, 'w') as temp:
			temp.write("> native " + str(self._start) + 'to' + str(self._end) + '\n')
			temp.write(self.seq_nuc() + '\n')
			temp.write("> pair " + str(OtherFrame._start) + 'to' + str(OtherFrame._end) + '\n')
			temp.write(OtherFrame.seq_nuc() + '\n')

		
		line = ['perl', '/mnt/gamma/user/sofya/tools/translatorx.pl', '-i', filename,\
		 '-o', filename[:-6], '-p', 'M', '-t', 'F', '-w', '1', '-c', '10']	
		process = subprocess.Popen(line)
		process.communicate()


		translator_aligned = reader.readfasta_tolist(filename[:-6] + ".nt_ali.fasta")
		

		# get a recursive list of file paths that matches pattern including sub directories
		fileList = glob.glob(filename[:-6] + "*")
		# Iterate over the list of filepaths & remove each file.
		for filePath in fileList:
		    try:
		        os.remove(filePath)
		    except OSError:
		        print("Error while deleting file")
		# os.remove(filename[:-6] + "*")


		# return pnps, [nd, N, sd, S] 
		dirname =  create_filename(flag = "paml_temp", extention = '', num = 15)
		os.mkdir('./' + dirname)
		os.system('cp ./codeml.clt ./' + dirname + '/codeml.clt')
		os.system("sed -i 's/.\/two_seqs.nuc/" + dirname + "\/two_seqs.nuc/g' " + dirname + '/codeml.clt')
		os.system("sed -i 's/codeml_2_seqs/" + dirname + "\/codeml_2_seqs/g' "+ dirname + '/codeml.clt')
		# seqfile = ./two_seqs.nuc * sequence data file name 
		# outfile = codeml_2_seqs * main result file)
 
		w, dN, N, dS, S= 100, 0, 0, 0, 0

		with open('/mnt/gamma/user/sofya/scripts/final/0pipeline/4_graphmark/' + dirname + '/two_seqs.nuc', 'w') as fout:
			fout.write('\t2\t' + str(len(translator_aligned[0])) + '\n')
			fout.write('num1' + '\n')
			fout.write(translator_aligned[0] + '\n')
			fout.write('num2' + '\n')
			fout.write(translator_aligned[1] + '\n')


		os.system('codeml /mnt/gamma/user/sofya/scripts/final/0pipeline/4_graphmark/' + dirname + '/codeml.clt  > ./' + dirname + '/codeml_2_seqs.out < codemlin.txt')

		with open(dirname + '/codeml_2_seqs', 'r') as fin:
			print('HEYss')
			for line in fin:
				if 'dN/dS=' in line: # t= 1.2947  S=   218.5  N=   582.5  dN/dS=  0.1324  dN = 0.1549  dS = 1.1693
					copy = line.replace('=', ' ')
					copy = copy.split()
					w = float(copy[7])
					dN = float(copy[9])
					dS = float(copy[11])
	
		with open(dirname + '/codeml_2_seqs.out', 'r') as fin:
			for line in fin:
				if 'dN*' in line: #             dN*=  0.15033 dS*=  1.27119 S* = 200.99 N* = 600.01
					copy = line.replace('=', ' ')
					copy = copy.split()
					S = float(copy[5])
					N = float(copy[7])
					break

		if (('paml_w' in self.out) == False):
			self.out['paml_w'] = list()
			self.out['paml_w_data'] = list()


		self.out['paml_w'].append(w)
		self.out['paml_w_data'].append([w, dN, N, dS, S])
		# pause = int(input())

		os.system('rm -r '+ dirname)
		
		return w, [dN, N, dS, S]




	# def biggest_ID(self):
	# 	identity_list = list()
	# 	if len(self._arrows) == 0:
	# 		return None

	# 	for arrow_to_frame in self._arrows:
	# 		s, u = self.find_pair_protein_id(arrow_to_frame)
	# 		if u > 5:
	# 			identity_list.append(float(s)/float(u))

	# 	a = max(identity_list)


	# 	self.out['biggest_prot_ID'] = max(identity_list)
	# 	# self.out['prot_IDs'] = ','.join(list(map(str,identity_list)))
	# 	return a

	def biggest_ungapped_ID(self):
		identity_list = list()
		identity_data = list()
		if len(self._arrows) == 0:
			self.out['biggest_prot_ID'] = None
			self.out['best_prot_ID'] = None
			self.out['all_id_data'] = None
			return None, None

		m_ungapped = 0
		m_same = 0 

		for arrow_to_frame in self._arrows:
			s, u = self.find_pair_protein_id(arrow_to_frame)
			print(s, u)
			if u > 5:
				identity_list.append(float(s)/float(u))
				identity_data.append(float(s)/float(u))
				identity_data.append(s)
				identity_data.append(u)
				if u > m_ungapped:
					m_ungapped = u 
					m_same = s
					self.reference_alignment = self.last_water_protein_alignment


		if identity_list == list():
				self.out['biggest_prot_ID'] = None
				self.out['best_prot_ID'] = None
				self.out['all_id_data'] = None
				return None, None

		self.out['biggest_prot_ID'] = max(identity_list)
		self.out['best_prot_ID'] = str(m_same/float(m_ungapped)) +',' + str(m_same) + ',' + str(m_ungapped)
		self.out['all_id_data'] =','.join(list(map(str,identity_data)))
		# self.out['prot_IDs'] = ','.join(list(map(str,identity_list)))
		return m_same, m_ungapped


	def fits_main_ID(self, reference_alignment, mainframe):
		# whatattt
		fits = False
		print(reference_alignment)
		
		if (self._start == mainframe._start) and (self._end == mainframe._end):
			fits = True
			self.out['fits_ID'] = fits
			self.out['bestfit_prot_ID'] ='main'
			return fits

		identity_list = list()
		best_fit = 0
		if self.out['all_id_data'] == None:
			self.out['fits_ID'] = fits
			self.out['bestfit_prot_ID'] = None
			return False
		idlist = list(map(float, self.out['all_id_data'].split(',')))

		for arrow_i in range(0, len(idlist), 3):
			my_meaning, same, ungapped = idlist[arrow_i:arrow_i + 3]


		# for arrow_to_frame in self._arrows:
		# 	same, ungapped = self.find_pair_protein_id(arrow_to_frame)

		# 	my_meaning = float(same)/float(ungapped)
		# 	identity_list.append(my_meaning)
			referencedata = cut_to_ids([reference_alignment[0], reference_alignment[1]], int(ungapped))
			# print(my_meaning)
			# print(referencedata)

			av = np.mean(np.array(referencedata))
			variance = np.var(np.array(referencedata))
			sigma = variance ** (0.5)

			try:
				sigmas_num = abs(av - my_meaning)/variance
				print(av)
				print(my_meaning)
				print(variance)
				print(sigmas_num)
			except ZeroDivisionError:
				if  (av - my_meaning == 0) :
					best_fit = my_meaning
					fits = True
				else: 
					continue 

			if self._start < mainframe._start:
				if ungapped < 10:
					continue
				elif ungapped > 60:
					if my_meaning > 0.4:
						best_fit = max(my_meaning, best_fit)
						fits = True
				elif sigmas_num <= before_sigmas[ungapped]:
					best_fit = max(my_meaning, best_fit)
					fits = True

			if self._end >  mainframe._end:
				if ungapped < 10:
					continue
				elif ungapped > 60:
					if my_meaning > 0.4:
						best_fit = max(my_meaning, best_fit)
						fits = True
				elif sigmas_num <=  after_sigmas[ungapped]:
					best_fit = max(my_meaning, best_fit)
					fits = True
		self.out['fits_ID'] = fits
		self.out['bestfit_prot_ID'] = best_fit
		# self.out['prot_IDs'] = ','.join(list(map(str, identity_list)))
		

		return fits




	
	def best_pnps(self): # a function for the main frame
		best_pnps = 100
		best_N = 0 
		best_counting_list = list()

		outline = str()

		for arrow_to_frame in self._arrows:
			pnps, counting_list = self.find_pair_pnps(arrow_to_frame)
			outline +=  str(pnps) + ',' + ','.join(list(map(str, counting_list))) + ',-,'
			nd2, N2, sd2, S2 = counting_list
			# if (pnps != None) and pnps < best_pnps:
			if (pnps != None) and best_N < N2:
				best_pnps = pnps
				best_N = N2
				best_counting_list = counting_list
		
		self.out['best_pnps'] = str(best_pnps)
		self.out['best_counting_list'] = ','.join(list(map(str, best_counting_list)))
		self.out['all_pnps'] = outline[:-1]
		return best_pnps, best_counting_list

	def best_paml_w(self): # a function for the main frame
		best_paml_w = 100
		best_N = 0 
		best_counting_list = list()

		outline = str()

		for arrow_to_frame in self._arrows:
			paml_w, counting_list = self.find_paml_w(arrow_to_frame)
			# print(paml_w)
			# pause = int(input())
			outline +=  str(paml_w) + ',' + ','.join(list(map(str, counting_list))) + ',-,'
			nd2, N2, sd2, S2 = counting_list
			# if ((paml_w != None) and (paml_w < best_paml_w)):
			if ((paml_w != None) and (best_N < N2)):
				best_paml_w = paml_w
				best_N = N2
				best_counting_list = counting_list
		# print(best_paml_w)
		# print(best_counting_list)
		# pause = int(input())
		self.out['best_paml_w'] = str(best_paml_w)
		self.out['best_paml_data'] = ','.join(list(map(str, best_counting_list)))
		self.out['all_paml_data'] = outline[:-1]
		
		return best_paml_w, best_counting_list


	def best_pnps_fit(self, nd1, N1, sd1, S1):
		pnps_fit_list = list()

		outline = str()

		for arrow_to_frame in self._arrows:
			pnps, counting_list = self.find_pair_pnps(arrow_to_frame)
			nd2, N2, sd2, S2 = counting_list
			if pnps != None:
				procent_of_second_bigger = get_p(nd1, N1, sd1, S1, nd2, N2, sd2, S2)
				outline += str(pnps) + ',' + ','.join(list(map(str, counting_list))) + ',' + str(procent_of_second_bigger)+ ',' 
				pnps_fit_list.append(procent_of_second_bigger)

		if len(pnps_fit_list) == 0:
			self.out['fit_best_pnps'] = None
			self.out['pnps_fit_list'] = None
			self.out['all_pnps'] = outline[:-1]
			return None
		
		self.out['fit_best_pnps'] = min(pnps_fit_list)
		self.out['pnps_fit_list'] = ','.join(list(map(str, pnps_fit_list)))
		self.out['all_pnps'] = outline[:-1]
		
		return min(pnps_fit_list)


	def best_paml_w_fit(self, nd1, N1, sd1, S1):
		paml_w_fit_list = list()

		outline = str()

		for arrow_to_frame in self._arrows:
			paml_w, counting_list = self.find_paml_w(arrow_to_frame)
			nd2, N2, sd2, S2 = counting_list
			if paml_w != None:
				procent_of_second_bigger = get_p(nd1, N1, sd1, S1, nd2, N2, sd2, S2)
				outline += str(paml_w) + ',' + ','.join(list(map(str, counting_list))) + ',' + str(procent_of_second_bigger)+ ',' 
				paml_w_fit_list.append(procent_of_second_bigger)

		if len(paml_w_fit_list) == 0:
			self.out['fit_best_paml_w'] = None
			self.out['paml_w_fit_list'] = None
			self.out['all_paml_data'] = outline[:-1]
			return None
		
		self.out['fit_best_paml_w'] = min(paml_w_fit_list)
		self.out['paml_w_fit_list'] = ','.join(list(map(str, paml_w_fit_list)))
		self.out['all_paml_data'] = outline[:-1]
		
		return min(paml_w_fit_list)




def MainFrame(Frame):
	def __init__(self, ):
		super(MainFrame, self).__init__()


def print_dict(frame_out): #idr _start _end length biggest_prot_ID prot_IDs fits_ID best_pnps best_counting_list fit_best_pnps pnps_fit_list all_pnps best_paml_w best_paml_data fit_best_paml_w paml_w_fit_list all_paml_data
	outlist = list()#idr _start _end length tomain is3 biggest_prot_ID best_prot_ID all_id_data fits_ID bestfit_prot_ID best_pnps best_counting_list fit_best_pnps pnps_fit_list all_pnps best_paml_w best_paml_data fit_best_paml_w paml_w_fit_list all_paml_data
	outlist.append(str(frame_out['_start']))
	outlist.append(str(frame_out['_end']))
	outlist.append(str(frame_out['length']))
	outlist.append(str(frame_out['tomain']))
	outlist.append(str(frame_out['is3']))
	try:
		outlist.append(str(frame_out['biggest_prot_ID']))
	except KeyError:
		outlist.append(str(None))
	try:
		outlist.append(str(frame_out['best_prot_ID']))
	except KeyError:
		outlist.append(str(None))
	try:
		outlist.append(str(frame_out['all_id_data']))
	except KeyError:
		outlist.append(str(None))
	try:
		outlist.append(str(frame_out['fits_ID']))
	except KeyError:
		outlist.append(str(None))
	try:
		outlist.append(str(frame_out['bestfit_prot_ID']))
	except KeyError:
		outlist.append(str(None))
	
	# try:
	# 	pri = str()
	# 	for pair in frame_out['prot_ids']:
	# 		pri += str(pair[0]) + ',' + str(pair[1]) + ','
	# 	outlist.append(str(pri[:-1]))
	# except KeyError:
	# 	outlist.append(str(None))
	# try:
	# 	outlist.append(str(frame_out['prot_IDs']))
	# except KeyError:
	# 	outlist.append(str(None))
	# try:
	# 	outlist.append(str(frame_out['fits_ID']))
	# except KeyError:
	# 	outlist.append(str(None))
	try:
		outlist.append(str(frame_out['best_pnps']))
	except KeyError:
		outlist.append(str(None))
	try:
		outlist.append(str(frame_out['best_counting_list']))
	except KeyError:
		outlist.append(str(None))
	try:
		outlist.append(str(frame_out['fit_best_pnps']))  
	except KeyError:
		outlist.append(str(None))
	try:
		outlist.append(str(frame_out['pnps_fit_list']))
	except KeyError:
		outlist.append(str(None))
	try:
		outlist.append(str(frame_out['all_pnps']))
	except KeyError:
		outlist.append(str(None))




	############################################
	# try:
	# 	outlist.append(str(frame_out['best_paml_w'])) 
	# except KeyError:
	# 	outlist.append(str(None))
	# try:
	# 	outlist.append(str(frame_out['best_paml_data']))
	# except KeyError:
	# 	outlist.append(str(None))
	# try:
	# 	outlist.append(str(frame_out['fit_best_paml_w']))  
	# except KeyError:
	# 	outlist.append(str(None))
	# try:
	# 	outlist.append(str(frame_out['paml_w_fit_list']))
	# except KeyError:
	# 	outlist.append(str(None))
	# try:
	# 	outlist.append(str(frame_out['all_paml_data']))
	# except KeyError:
	# 	outlist.append(str(None))
	############################################




	return '\t'.join(outlist)




"""making a list of frames - from marking module steal 
the stopstop list maker """

stops = ['TAA', 'TAG']# , 'TGA']

def create_frame_list(seq): #and we need the first and the last. AND we dont need to differ them by frame
	stops = ['TAA', 'TAG']
	Frs = list()
	
	i = 0
	for j in range(0, len(seq), 3):
		if (seq[j:j+3] in stops) or (j  >= len(seq) - 3):
			fr = [i, j]
			i = j + 3
			Frs.append(fr)

	i = 1
	for j in range(1, len(seq), 3):
		if (seq[j:j+3] in stops) or (j  >= len(seq) - 3):
			fr = [i, j]
			i = j + 3
			Frs.append(fr)

	i = 2
	for j in range(2, len(seq), 3):
		if (seq[j:j+3] in stops) or (j  >= len(seq) - 3):
			fr = [i, j]
			i = j + 3
			Frs.append(fr)

	return Frs


#THIS sequence is to be cleaned if no exceptions are ever raised
def get_inalignmment(alt_seq, posstart, posend):
	# print(posstart, posend)
	alstart = 0 
	alend = 0 
	filled = 0 
	now = -1 
	
	while ((filled <= posstart) and (now < len(alt_seq))) :
		now += 1
		if now == len(alt_seq):
			break
		try: 
			if alt_seq[now] != '-':
				filled += 1
		except IndexError:
			print(now)
			print(len(alt_seq))
			print(posend)
			print(posstart)

			pause = int(input())
	
	alstart = now
	# print(alstart)

	while ((filled <= posend) and  (now < len(alt_seq))) :
		now += 1
		if now == len(alt_seq):
			break
		try: 
			if alt_seq[now] != '-':
				filled += 1
		except IndexError:
			print('ERROR IN GET_INALIGNMENT')
			print(now)
			print(len(alt_seq))
			print(posend)
			print(posstart)
			pause = int(input())

	alend = now
	# print(alstart, alend)
	return alstart, alend



def fr_length(fr):
	return (fr[1] - fr[0])//3

def filter_by_threshold(threshold, Frs):

	Nice = list()

	for fr in Frs: 
		if fr_length(fr) >= threshold:
			Nice.append(fr)

	return Nice

def turn_Frs_to_alignmentFrs(Frs, alt_seq):
	AlFrs = list()
	for fr in Frs:
		s, e = get_inalignmment(alt_seq, fr[0], fr[1])
		AlFrs.append([s, e])
	return AlFrs

def make_object_Frs(AlFrs, alt_seq):
	objFrs = list()
	for fr in AlFrs:
		objFrs.append(Frame(start = fr[0], end = fr[1], parent = alt_seq))
	return objFrs


def byStart(elem):
    return elem._start
def byEnd(elem):
    return elem._end



def intercept(fr1, fr2): #IN BP
	s1, e1 = fr1._start, fr1._end
	s2, e2 = fr2._start, fr2._end
	if (e1 <= s2) or (e2 <= s1):
		return 0
	s = max(s1, s2)
	e = min(e1, e2)
	if (s > e):
		return 0 
	seq1 = fr1._parent
	seq2 = fr2._parent
	same, ungapped = find_id(seq1[s:e], seq2[s:e])
	return ungapped



def Build_Arrows(Frs_from, Frs_to, threshold_ungapped):
	start_sort = sorted(Frs_to, key=lambda x: x._start)
	# end_sort = sorted(Frs_to, key=bySecond)

	for fr1 in Frs_from:
		fr1._arrows = list()
		for fr2 in start_sort:
			if fr2._start > fr1._end:
				break
			if fr2._end < fr1._start:
				continue
			if intercept(fr1, fr2)//3 > threshold_ungapped:
				fr1.add_arrow(fr2)


class Transcript(object):
	def __init__(self, idr, bp_align = None):
		self._idr = idr
		self._species = idr[:4]
		self.bp_align = bp_align

		# if bp_align == None: 
		# 	self.get_alignment_fromfile()

		



	# def get_alignment_fromfile(self, directory='/mnt/gamma/user/sofya/scripts/final/homo_bp_rev_alignwater/*'): 
	# 	args = ['gawk', '-e ', '"/' + str(self._idr) + '/{getline; print; exit}"', directory, '>', 'temp.txt']
	# 	cmd = ' '.join(args)

	# 	os.system(cmd)
	# 	with open('temp.txt', 'r') as fin:

	# 		e = fin.readline()

	# 	if e == str():
	# 		print('EXCEPTION from Transcript.get_sequence_fromfile:\
	# 			\nNo sequence with {0} identificator in {1} directory'.format(self._idr, directory))
	# 		return 
	# 	self.bp_align = e[:-1]



	def frames_preprocess(self, threshold): # Here is diffenent from NEW!
		# global filtstats1
		self.frames = create_frame_list(clean(self.bp_align))
		# filtstats1.write(str(len(self.frames)) + '\t')
		self.frames = filter_by_threshold(threshold, self.frames)
		# filtstats1.write(str(len(self.frames)) + '\t')
		self.frames = turn_Frs_to_alignmentFrs(self.frames, self.bp_align)
		
		self.frames = make_object_Frs(self.frames, self.bp_align)
		


	def filter_frames_p1(self, Ortholog, threshold_ungapped):
		Build_Arrows(self.frames, Ortholog.frames, threshold_ungapped)
		print(len(self.frames ))
		# pause  = int(input())
		
		longest_ungapped = 0

		for fr in self.frames:
			e = fr.get_length()
			s, u = fr.biggest_ungapped_ID()
			print(s, u)
			if u == None:
				continue
			
			if u > longest_ungapped:
				longest_ungapped = u
				self.mainframe = fr

		outside_of_main = list()

		for fr in self.frames:
			s = fr._start
			e = fr._end
			# print()
			if (s > self.mainframe._start) and (e < self.mainframe._end): # if the frame is fully inside the mainframe - delete it. Though we dont filter out the minframe itself here!!!!
				continue
			elif (s < self.mainframe._end) and (e > self.mainframe._end):
				fr.out['tomain'] = 3
				outside_of_main.append(fr)

				add_s = self.mainframe._end + 1
				while len(clean(self.bp_align[add_s:e]))%3 != 0:
					add_s += 1

				also = Frame(start = add_s, end = e, parent = self.bp_align, tomain = 3 , is3 = True)
				l = also.get_length()
				s, u = also.biggest_ungapped_ID()
				outside_of_main.append(also)
			elif (s == self.mainframe._start) and (e == self.mainframe._end):
				fr.out['tomain'] = 'main'
				outside_of_main.append(fr)
			else:
				fr.out['tomain'] = 5
				outside_of_main.append(fr)

		if (self.mainframe in outside_of_main) == False:
			print('WARNING::transcript class::frames_preprocess:: something went wrong I threw out the main frame please rewrite the for cycle ')
			pause = int(input())

		self.frames = outside_of_main
		# print(len(self.frames ))
		# pause  = int(input())



	def filter_frames_p2(self, Ortholog, threshold_ungapped):
		Build_Arrows(self.frames, Ortholog.frames, threshold_ungapped)
		# print(len(self.frames ))
		# pause  = int(input())
		# filtstats1.write(str(len(self.frames)) + '\n')



	# def filter_frames(self, Ortholog, threshold_ungapped):
		# global filtstats2
		global logfile
		# Build_Arrows(self.frames, Ortholog.frames, threshold_ungapped)
		forstat = 0
		left = list()
		# self.mainframe
		# a, b = self.mainframe.find_pair_protein_id(Ortholog.mainframe)
		# self.mainframe.biggest_ID()
		# self.mainframe.out['biggest_ID'] = str(float(a)/float(b))
		reference_alignment = self.mainframe.reference_alignment 
		
		pnps, counting_list  = self.mainframe.best_pnps()
		

		nd1, N1, sd1, S1 = counting_list


		############################################
		# w, paml_counting_list  = self.mainframe.best_paml_w()
		# nd1p, N1p, sd1p, S1p = counting_list
		############################################
		

		# ID_threshold = 0.7

		if len(self.frames) == 0:
			return None



		# outside_of_main = list()

		for fr in self.frames:
			s = fr._start
			e = fr._end
			if (s > self.mainframe._start) and (e < self.mainframe._end): # if the frame is fully inside the mainframe - delete it. Though we dont filter out the minframe itself here!!!!
				fr.out['tomain'] = 'inside'
				# continue
			elif (s == self.mainframe._start) and (e == self.mainframe._end):
				fr.out['tomain'] = 'main'
			# else:
			# 	fr.out['tomain'] = 'outside'
				# outside_of_main.append(fr)


		for fr in self.frames:
			fits_ID = fr.fits_main_ID(reference_alignment, self.mainframe) #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			
			best_pnps_fit = fr.best_pnps_fit(nd1, N1, sd1, S1)
			
			############################################
			# best_paml_w_fit = fr.best_paml_w_fit(nd1p, N1p, sd1p, S1p)
			############################################
			
			logfile.write(self._idr + '\t')
			logfile.write(print_dict(fr.out))
			logfile.write('\n')

			if (fr.out['length']//3 > 60):
				left.append(fr)
			elif (fits_ID == True):
				# forstat += 1
				if (best_pnps_fit != None) and (best_pnps_fit < 0.8):
					left.append(fr)

		# filtstats2.write(str(forstat) + '\t' + str(len(left)) + '\n')

		return left



	def get_coding(self, Ortholog, intercept_treshold): #Ortholog must be preprocessed
		left = self.filter_frames_p1(Ortholog, intercept_treshold)
		# ort = Ortholog.filter_frames_p1(self, intercept_treshold)
		left = self.filter_frames_p2(Ortholog, intercept_treshold)
		left = sorted(left, key=lambda x: x._start)

		# if 


		# with open('log5.txt', 'a') as fout:	
		# 	for fr in left:
		# 		fout.write(self._idr +  '\t' + fr.output() + '\n')

		if left == None:
			return '[]'
		out = str()

		for fr in left:
			out += '[' + str(fr._start) + ', ' + str(fr._end) + '] '
		
		return out









"""
HERE STARTS THE MAIN PART
here we sort by length and leave only those that have 
a length more than something"""

"""than we shall make the arrow connections from one 
sequence to tthe other"""
"""for each of the sequences we get the BEST protein 
ID and filter those that have a low one """
"""than the best pnps fit and sort them """
"""than we are left with the best frames """
"""and try to connect them"""




# directory = '/mnt/gamma/user/sofya/scripts/final/0pipeline/12_golden/nuc_gold/'
directory = '/mnt/gamma/user/sofya/scripts/final/0pipeline/nuc_alignments/'
files  = os.listdir(directory)

# # # files = ['/mnt/gamma/user/sofya/scripts/final/homo_bp/hom_cras_c15297_g1_i1.fasta']

parser = OptionParser()
parser.add_option("-f", "--flag", help="for output files")
parser.add_option("-s", "--interval_start", help="some interval of files to look at")
parser.add_option("-e", "--interval_end", help="some interval of files to look at") #1469 files all in all 
opt, args = parser.parse_args()


flag =  opt.flag
parts = int(opt.interval_start)
parte = int(opt.interval_end)




# files = ['hom_cras_c2265_g1_i1.fasta']
# flag = 'test2'
# files = ['hom_cras_c2265_g1_i1.fasta',
# 'hom_cras_c25296_g1_i1.fasta',
# 'hom_cras_c19750_g1_i1.fasta',
# 'hom_cras_c13999_g1_i1.fasta',
# 'hom_cras_c25771_g1_i1.fasta',
# 'hom_cras_c5734_g1_i1.fasta',
# 'hom_cras_c2551_g1_i1.fasta',
# 'hom_cras_c5314_g1_i1.fasta',
# 'hom_cras_c3177_g1_i1.fasta',
# 'hom_cras_c8927_g1_i1.fasta',
# 'hom_cras_c7010_g1_i1.fasta',
# 'hom_cras_c24489_g1_i1.fasta',
# 'hom_cras_c7233_g1_i1.fasta',
# 'hom_cras_c5431_g1_i1.fasta',
# 'hom_cras_c4608_g1_i1.fasta',
# 'hom_cras_c6045_g1_i1.fasta',
# 'hom_cras_c24279_g1_i1.fasta',
# 'hom_cras_c6491_g1_i1.fasta',
# 'hom_cras_c5570_g1_i1.fasta',
# 'hom_cras_c29525_g1_i1.fasta',
# 'hom_cras_c8091_g1_i1.fasta',
# 'hom_cras_c17093_g1_i1.fasta',
# 'hom_cras_c30509_g1_i1.fasta',
# 'hom_cras_c1282_g1_i1.fasta',
# 'hom_cras_c7093_g1_i1.fasta',
# 'hom_cras_c16907_g1_i1.fasta']



def filter_paralogues(bp_dict):
	animals = list()
	para = list()
	new_dict = dict()

	for key in bp_dict.keys():
		animal  = key[:4]
		if animal in animals:
			para.append(animal)
		animals.append(key[:4])

	for key in bp_dict:
		animal  = key[:4]
		if animal in para:
			continue
		else:
			new_dict[key] = bp_dict[key]

	return new_dict



logfile = open('./logs_from_graphmarker' + flag + '.txt', 'a') #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# files = ['hom_cras_c7093_g1_i1.fasta',
# 'hom_raik_c8008_g1_i1.fasta',
# 'hom_cras_c1282_g1_i1.fasta',
# 'hom_cras_c2551_g1_i1.fasta',
# 'hom_cras_c25771_g1_i1.fasta',
# 'hom_eury_eury_c12004_g1_i1.fasta',
# 'hom_foca_c10387_g2_i1.fasta',
# 'hom_eury_eury_c15759_g1_i1.fasta',
# 'hom_cras_c19750_g1_i1.fasta',
# 'hom_cras_c24489_g1_i1.fasta']



for filename in files[parts:parte]:
# for filename in files:


	path = directory + filename
	bp_dict = reader.readfasta_todictionary(path)
	bp_dict = filter_paralogues(bp_dict)
	idrs = bp_dict.keys()

	transcript_dict = dict()
	longest = str()
	llen = 0 
	secondlongest = str()
	sllen = 0 

	for idr in bp_dict:
		transcript_dict[idr] = Transcript(idr, bp_align =bp_dict[idr])
		transcript_dict[idr].frames_preprocess(10)
		if len(clean(bp_dict[idr])) > llen:
			secondlongest = longest
			sllen = llen
			longest = idr
			llen = len(clean(bp_dict[idr]))
		elif len(clean(bp_dict[idr])) > sllen:
			secondlongest = idr
			sllen = len(clean(bp_dict[idr]))

	for idr in idrs:

		if idr != longest:
			ortholog =  transcript_dict[longest]
		else:
			ortholog = transcript_dict[secondlongest]

		coding = transcript_dict[idr].get_coding(ortholog, 5)

logfile.close()
