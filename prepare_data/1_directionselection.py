# sofya.gaydukova@gmail.com 02/03/2021
#
# Select a direction of a sequence. Use dN/dS ratio. 
# Choose the longest frame in the + and - stranis and compair their dN/dS.
#
#
#


import re
import subprocess
import shlex
import os
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader
import marking_module_v3 as mark
import bootstrap_tables_v4 as boot
import ds
# import dndsframecutter as cut
from optparse import OptionParser



parser = OptionParser()
# parser.add_option("-i", "--identificators", help="commseqids", default='foca_c7473_g1_i1')
parser.add_option("-b", "--front5", help="threshold for the length of the 5'end pieces", default=1000)
parser.add_option("-s", "--startorf", help="threshold for the length of the first detected orf", default=30)
parser.add_option("-e", "--end3", help="threshold for the length of the 3'end pieces", default=1000)


parser.add_option("-d", "--directory", help="directory with clean orthogoups")


opt, args = parser.parse_args()








# idrs = str(opt.identificators)
# idrs = list(idrs.split(','))

front = int(opt.front5)
end = int(opt.end3)
start = int(opt.startorf)




# from optparse import OptionParser

# parser = OptionParser()
# parser.add_option("-d", "--directory", help="directory with clean orthogoups")
# opt, args = parser.parse_args()


directory =  opt.directory

# directory = '/mnt/gamma/user/sofya/scripts/final/homo_bp_rev_alignedwater/'
files  = os.listdir(directory)



def get_inalignmment(ald_seq, posstart, posend):
	# print(posstart, posend)
	alstart = 0 
	alend = 0 
	filled = 0 
	now = -1 
	
	while ((filled <= posstart) and (now < len(ald_seq))) :
		now += 1
		if now == len(ald_seq):
			break
		try: 
			if ald_seq[now] != '-':
				filled += 1
		except IndexError:
			print(now)
			print(len(allseq))
			print(posend)
			print(posstart)

			pause = int(input())
	
	alstart = now
	# print(alstart)

	while ((filled <= posend) and  (now < len(ald_seq))) :
		now += 1
		if now == len(ald_seq):
			break
		try: 
			if ald_seq[now] != '-':
				filled += 1
		except IndexError:
			print('ERROR IN GET_INALIGNMENT')
			print(now)
			print(len(allseq))
			print(posend)
			print(posstart)
			pause = int(input())

	alend = now
	# print(alstart, alend)
	return alstart, alend




class Transcript(object):
	def __init__(self, idr):
		self._idr = idr
		self._species = idr[:4]
		self.codings = dict()
		self.saved = dict()

	# need someone to look at please 
	def get_sequence_fromfile(self, directory='./orthogroups/*'): 
		print("FOR NORMAL")
		idr_for_gawk = self._idr.replace('+', '\\+')
		if directory[-1] != '*':
			args = ['gawk', '-e ', '"/' + str(idr_for_gawk ) + '/{getline; print; exit}"', directory + '*', '>', 'temp.txt']
		else:
			args = ['gawk', '-e ', '"/' + str(idr_for_gawk ) + '/{getline; print; exit}"', directory, '>', 'temp.txt']
		# args = ['gawk', '-e ', '"/' + str(idr) + '/{getline; print; exit}"', './hello/sdf.fasta']
		# args = ['awk', '-e ', '"/' + str(idr) ./orthogroups/*.fasta']
		cmd = ' '.join(args)

		os.system(cmd)
		print(cmd)
		with open('temp.txt', 'r') as fin:
			e = fin.readline()


		if e == str():
			print('EXCEPTION from Transcript.get_sequence_fromfile:\
				\nNo sequence with {0} identificator in {1} directory'.format(self._idr, directory))
			return 
		self.bp_seq = e[:-1]

	def get_complement(self):
		return mark.complementary(self.bp_seq)
	
	def get_coding(self, front = 1000, end=1000, start=40):
		# print(self.bp_seq)
		cod, s = mark.get_coding_only_straight(self.bp_seq, front, end, start)
		# print(cod,s )
		return cod

	def framecut(self, front = 1000, end=1000, start=40):
		breakinfo = list()

		coding = self.get_coding(front, end, start)
		cut = str()
			
		if len(coding) == 0:
			breakinfo  = 'not_coding'
			return (breakinfo, str())

		stopplace = 0 

		for b in range(0, len(coding), 2):
			cut +=  self.bp_seq[coding[b]:coding[b+1]]
			stopplace += coding[b + 1] - coding[b]
			breakinfo.append(stopplace)
		self.codings[','.join(list(map(str, coding)))] = tuple([breakinfo, cut])
		return tuple([breakinfo, cut])

	def add_to_file_framecut(self, filepath, mode='a', front = 1000, end=1000, start=40):
		breakinfo, cut = self.framecut(front , end, start)
		with open(filepath, mode) as fout:
			fout.write('> ' + self._idr + ' ' + ','.join(list(map(str, breakinfo))) + '\n')
			fout.write(cut + '\n')

	def saved_string(self):
		st = str()
		st = self._idr
		for elem in ['coding', 'breakinfo', 'id_water_cut', 'id_water_uncut', 'id_codon_cut', 'pnps_codon_cut']:
			try:
				st += '\t' + ','.join(list(map(str, self.saved[elem])))
			except TypeError:
				st += '\t' + str(self.saved[elem])
		print(st)
		return st 

class RevTranscript(Transcript):
	def __init__(self, idr):
		# super(RevTranscript, self).__init__(idr)
		super(RevTranscript, self).__init__(idr)


	def get_sequence_fromfile(self, directory='./orthogroups/*'): 
		print("FOR REVERSED")
		idr_for_gawk = self._idr.replace('+', '\\+')
		if directory[-1] != '*':
			args = ['gawk', '-e ', '"/' + str(idr_for_gawk ) + '/{getline; print; exit}"', directory + '*', '>', 'temp.txt']
		else:
			args = ['gawk', '-e ', '"/' + str(idr_for_gawk ) + '/{getline; print; exit}"', directory, '>', 'temp.txt']

		cmd = ' '.join(args)
		os.system(cmd)
		print(cmd)
		with open('temp.txt', 'r') as fin:

			e = fin.readline()
		if e == str():
			print('EXCEPTION from Transcript.get_sequence_fromfile:\
				\nNo sequence with {0} identificator in {1} directory'.format(self._idr, directory))
			return 
		self.bp_seq = mark.complementary(e[:-1])


class OrthoPair(object):
	def __init__(self, idr_list=None, directory = './orthogroups/*'):
		self._idr_list = idr_list

		self._transcript_list = list([Transcript(self._idr_list[0]), Transcript(self._idr_list[1])])
		self._species = [idr_list[0][:4], idr_list[1][:4]]
		self._bp_seqs = list()
		for tr in self._transcript_list:
			tr.get_sequence_fromfile(directory)
			self._bp_seqs.append(tr.bp_seq)
			print("ss", tr._idr, tr.bp_seq[:10] , tr.bp_seq[-10:])
		self._bp_dict = {self._idr_list[0]:self._bp_seqs[0], self._idr_list[1]:self._bp_seqs[1]}
		# self._bp_dict
		# print(self._bp_dict)


	def pairwise_codon_alignment(self, front = 1000, end = 1000, start = 40):

	# 	# outline1 = k[0] + '\t' +  k[1] + '\t' + str(len(bp_dict[k[0]]))  + '\t'+ str(len(bp_dict[k[0]]))
		for tr in self._transcript_list:
			tr.saved['fullen'] = len(tr.bp_seq)
			tr.saved['coding'] = tr.get_coding()



		cut_dict = dict()
		info_dict = dict()

		for tr in self._transcript_list:
			breakinfo, cut = tr.framecut(front, end, start)
			cut_dict[tr._idr] = cut
			tr.saved['breakinfo'] = breakinfo
			tr.saved['cut'] = cut
			info_dict[tr._idr] = breakinfo
			if breakinfo == 'not_coding':
				print('ERROR in pairwise_codon_alignment - oni is not coding!!')
				return 'not_coding', 'not_coding', 'not_coding'

		water_aligned_cut = boot.getanalignment('forpairwise.fasta', cut_dict)
		same, length = boot.find_id(water_aligned_cut[self._idr_list[0]], water_aligned_cut[self._idr_list[1]])
		identity = float(same) / float(length)
		for tr in self._transcript_list:
			tr.saved['id_water_cut'] = identity
			tr.saved['water_cut'] = water_aligned_cut[tr._idr]

		water_aligned_uncut = boot.getanalignment('forpairwise.fasta', self._bp_dict)
		same, length = boot.find_id(water_aligned_uncut[self._idr_list[0]], water_aligned_uncut[self._idr_list[1]])
		identity = float(same) / float(length)
		for tr in self._transcript_list:
			tr.saved['id_water_uncut'] = identity
			tr.saved['water_uncut'] = water_aligned_uncut[tr._idr]
		# print('---------------------===============================')
		# print(identity)

		# with open('/mnt/gamma/user/sofya/scripts/final/dnds/2_coding_parts_unaligned_forone.fasta', 'w') as temp:
		with open('temp.fasta', 'w') as temp:
			for tr in self._transcript_list:
				temp.write('> ' + tr._idr + '\n')
				temp.write(tr.saved['cut'] + '\n')
		
		line = ['perl', '/mnt/gamma/user/sofya/tools/translatorx.pl', '-i', \
		'temp.fasta',\
		 '-o', 'temp2', \
		 '-p', 'M', '-t', 'F', '-w', '1', '-c', '10']	
		process = subprocess.Popen(line)
		process.communicate()

		translator_aligned_cut = reader.readfasta_tolist('temp2.nt_ali.fasta')
		translator_aligned_cut_dict = reader.readfasta_todictionary('temp2.nt_ali.fasta')
		
		nd, N, sd, S = ds.getdnds(translator_aligned_cut[0], translator_aligned_cut[1])
		try: 
			pnps = float(nd)/ float(N) / float(sd) * float(S)

		except ZeroDivisionError:
			pnps = 'NULL'
		output = pnps

		same, length = boot.find_id(translator_aligned_cut[0], translator_aligned_cut[1])
		identity = float(same) / float(length)
		
		for tr in self._transcript_list:
			tr.saved['id_codon_cut'] = (identity, same, length)
			tr.saved['codon_cut'] = translator_aligned_cut_dict['_' + tr._idr]
			tr.saved['pnps_codon_cut'] = ( nd, N, sd, S, pnps)
 
		

		outputline1 = self._transcript_list[0].saved_string()
		outputline2 = self._transcript_list[1].saved_string()
		print(outputline1)
		print(outputline2)

		# for tr in self._transcript_list:
		# 	tr.saved['pieces'] = str()
		# 	for index in range(len(tr.saved['breakinfo'])):
		# 		somestop = tr.saved['breakinfo'][index]
		# 		if index == 0:
		# 			before = 0
		# 		else:
		# 			before = tr.saved['breakinfo'][index - 1]
		# 		piece_len = somestop - before

		# 		al_before, al_somestop = get_inalignmment(translator_aligned_cut_dict['_' +tr._idr], before, somestop)

		# 		sd, S, nd, N = ds.getdnds(translator_aligned_cut[0][al_before:al_somestop], translator_aligned_cut[1][al_before:al_somestop])
		# 		same, length = boot.find_id(translator_aligned_cut[0][al_before:al_somestop], translator_aligned_cut[1][al_before:al_somestop])			
		# 		try: 
		# 			pnps = float(nd)/ float(N) / float(sd) * float(S)
		# 		except ZeroDivisionError:
		# 			pnps = 'NULL'

		# 		tr.saved['pieces'] += '_' + ','.join(list(map(str, [piece_len, sd, S, nd, N, pnps, same, length, identity])))
		return output, outputline1, outputline2

class RevOrthoPair(OrthoPair):
	def __init__(self, idr_list, directory = './orthogroups/*'):

		self._idr_list = idr_list
		self._transcript_list = list([RevTranscript(self._idr_list[0]), RevTranscript(self._idr_list[1])])


		self._species = [idr_list[0][:4], idr_list[1][:4]]
		self._bp_seqs = list()
		for tr in self._transcript_list:
			tr.get_sequence_fromfile()
			self._bp_seqs.append(tr.bp_seq)
			print("rr", tr._idr, tr.bp_seq[:10] , tr.bp_seq[-10:])

		self._bp_dict = {self._idr_list[0]:self._bp_seqs[0], self._idr_list[1]:self._bp_seqs[1]}
		# print(self._bp_dict)

class FirstRevOrthoPair(OrthoPair):
	def __init__(self, idr_list, directory = './orthogroups/*'):

		self._idr_list = idr_list
		self._transcript_list = list([RevTranscript(self._idr_list[0]), Transcript(self._idr_list[1])])


		self._species = [idr_list[0][:4], idr_list[1][:4]]
		self._bp_seqs = list()
		for tr in self._transcript_list:
			tr.get_sequence_fromfile()
			self._bp_seqs.append(tr.bp_seq)
			print("rs", tr._idr, tr.bp_seq[:10] , tr.bp_seq[-10:])
		self._bp_dict = {self._idr_list[0]:self._bp_seqs[0], self._idr_list[1]:self._bp_seqs[1]}

class SecondRevOrthoPair(OrthoPair):
	def __init__(self, idr_list, directory = './orthogroups/*'):

		self._idr_list = idr_list
		self._transcript_list = list([Transcript(self._idr_list[0]), RevTranscript(self._idr_list[1])])


		self._species = [idr_list[0][:4], idr_list[1][:4]]
		self._bp_seqs = list()
		for tr in self._transcript_list:
			tr.get_sequence_fromfile()
			self._bp_seqs.append(tr.bp_seq)
			print("sr", tr._idr, tr.bp_seq[:10] , tr.bp_seq[-10:])
		self._bp_dict = {self._idr_list[0]:self._bp_seqs[0], self._idr_list[1]:self._bp_seqs[1]}



# ids = ['harp_c33973_g2_i1', 'octo_c45705_g1_i1']

# def output(tr):
# 	print('seq',tr.bp_seq)
# 	print('fullen', tr.saved['fullen'])
# 	print('water_uncut', tr.saved['water_uncut'])
# 	print('id_water_uncut', tr.saved['id_water_uncut'])

# 	print('coding', tr.saved['coding'])
# 	print('breakinfo', tr.saved['breakinfo'])
# 	print('cut', tr.saved['cut'])
# 	print('water_cut', tr.saved['water_cut'])
# 	print('id_water_cut', tr.saved['id_water_cut'])

# 	print('codon_cut', tr.saved['codon_cut'])
# 	print('id_codon_cut', tr.saved['id_codon_cut'])
# 	print('pnps_codon_cut', tr.saved['pnps_codon_cut'])
# 	print('pieces', tr.saved['pieces'])
# 	print('\n')


# ids = ['cras_c6591_g1_i1', 'foca_c7473_g1_i1']

# # al  = OrthoPair(ids)
# # print(al._bp_seqs)  

# st_st1 = OrthoPair(ids)
# st_st1.pairwise_codon_alignment()

# for tr in st_st1._transcript_list:
# 	output(tr)
	# print(tr.saved)



# st_st2 = RevOrthoPair(ids)
# st_st2.pairwise_codon_alignment()

# for tr in st_st2._transcript_list:
# 	output(tr)
# 	# print(tr.saved)



# st_st3 = HalfRevOrthoPair(ids)
# st_st3.pairwise_codon_alignment()

# for tr in st_st3._transcript_list:
# 	output(tr)
# 	# print(tr.saved)



# st_st4 = HalfRevOrthoPair([ids[1], ids[0]])
# st_st4.pairwise_codon_alignment()

# for tr in st_st4._transcript_list:
# 	output(tr)
	# print(tr.saved)





# files = ['/mnt/gamma/user/sofya/scripts/final/homo_bp/hom_cras_c15297_g1_i1.fasta']



def choose_direction(ref, idr2, directory):
	st_st = OrthoPair([ref, idr2], directory)
	ss, ss1, ss2 = st_st.pairwise_codon_alignment()

	re_st = FirstRevOrthoPair([ref, idr2], directory)
	rs, rs1, rs2 = re_st.pairwise_codon_alignment()

	st_re = SecondRevOrthoPair([ref, idr2], directory)
	sr, sr1, sr2 = st_re.pairwise_codon_alignment()

	re_re = RevOrthoPair([ref, idr2], directory)
	rr, rr1, rr2 = re_re.pairwise_codon_alignment()

	if ss == min(ss, rs, sr, rr):
		return "ss", idr2 + '--' + str(ss) + '--' + ss2 + '--' + rs2  + '--' + sr2  + '--' + rr2, ref + '--' + str(ss) + '--' + ss1 + '--' + sr1  + '--' + rs1  + '--' + rr1
	elif rs == min(ss, rs, sr, rr):
		return "rs", idr2 + '--' + str(rs) + '--' + rs2 + '--' + ss2  + '--' + rr2  + '--' + sr2, ref + '--' + str(rs) + '--' + rs1 + '--' + rr1  + '--' + ss1  + '--' + sr1
	elif sr == min(ss, rs, sr, rr):
		return "sr", idr2 + '--' + str(sr) + '--' + sr2 + '--' + rr2  + '--' + ss2  + '--' + rs2, ref + '--' + str(sr) + '--' + sr1 + '--' + ss1  + '--' + rr1  + '--' + rs1
	elif rr == min(ss, rs, sr, rr):
		return "rr", idr2 + '--' + str(rr) + '--' + rr2 + '--' + sr2  + '--' + rs2  + '--' + ss2, ref + '--' + str(rr) + '--' + rr1 + '--' + rs1  + '--' + sr1  + '--' + ss1







already = os.listdir('./ortho_turned/')
print(already)


info = open('infofile_2.txt', 'w')


# already = os.listdir('/mnt/gamma./orthogroups/')

for filename in files:
	# if filename in already:
	# 	continue
	path = directory + filename
	bp_dict = reader.readfasta_todictionary(path)
	bp_dict = reader.filter_paralogues(bp_dict)

	ref = reader.choose_reference(bp_dict)
	tag = 0 
	outfilepath = './ortho_turned/' + filename
	if filename in already:
		print('yes')
		continue
	outfile = open(outfilepath, 'w')
	for idr in bp_dict:
		if idr == ref:
			continue
		dire, idr_st, ref_st = choose_direction(ref, idr, directory)
		info.write(idr_st + '\n')
		if tag == 0:
			info.write(ref_st + '\n')

			if dire[0] == 's' :
				outfile.write("> " + ref + '\n')
				outfile.write(bp_dict[ref] + '\n')
			elif dire[0] == 'r' :
				outfile.write("> " + ref + '\n')
				outfile.write(mark.complementary(bp_dict[ref]) + '\n')
			tag = 1 
		
		if dire in ["ss","rs"] :
			outfile.write("> " + idr + '\n')
			outfile.write(bp_dict[idr] + '\n')
		elif dire in ["sr","rr"] :
			outfile.write("> " + idr + '\n')
			outfile.write(mark.complementary(bp_dict[idr]) + '\n')

		# pause = int(input())

	outfile.close()






# seqs = list(getseq(idr[0]), getseq(idr[1]))
# bp_dict = {idr[0]:seqs[0], idr[1]:seqs[1]}



# # nothing here reverses anything!
# def pairwise_codon_alignment(bp_dict, front, end, start):



# 	# if 'TTTTTTTT' == seqs[0][:8]:
# 	# 	to_write[k[0]] = mark.complementary(bp_dict[k[0]])
# 	# if 'TTTTTTTT' == bp_dict[k[1]][:8]:
# 	# 	to_write[k[1]] = mark.complementary(bp_dict[k[1]])
# 	k = bp_dict.keys()

# 	# outline1 = k[0] + '\t' +  k[1] + '\t' + str(len(bp_dict[k[0]]))  + '\t'+ str(len(bp_dict[k[0]]))
# 	outline1 = k[0] + '\t' + str(len(bp_dict[k[0]]))
# 	outline2 = k[1] + '\t' + str(len(bp_dict[k[1]]))

# 	end = 5
# 	for end in range(5, 6):


# #------------------------------------------ align ----------------------------------------------------------
# 		infodict = cut.framecut_only_straight(to_write, '/mnt/gamma/user/sofya/scripts/final/dnds/2_coding_parts_unaligned_endscopy.fasta'\
# 			, front, end, start)
		
# 		if infodict[k[0]] == 'not_coding' or infodict[k[1]] == 'not_coding':
# 			outline1 += '\tnc\tnc\tnc\tnc\tnc\tnc\tnc\tnc\tnc'
# 			outline2 += '\tnc\tnc\tnc\tnc\tnc\tnc\tnc\tnc\tnc'
# 			continue 


# 		line = ['perl', '/mnt/gamma/user/sofya/tools/translatorx.pl', '-i', \
# 		'/mnt/gamma/user/sofya/scripts/final/dnds/2_coding_parts_unaligned_endscopy.fasta',\
# 		 '-o', '/mnt/gamma/user/sofya/scripts/final/dnds/3_aligned_seqs_endscopy', \
# 		 '-p', 'M', '-t', 'F', '-w', '1', '-c', '10']	
# 		process = subprocess.Popen(line)
# 		process.communicate()


# 		translated = reader.readfasta_tolist('/mnt/gamma/user/sofya/scripts/final/dnds/3_aligned_seqs_endscopy.nt_ali.fasta')
# 		sd, S, nd, N = ds.getdnds(translated[0], translated[1])


# 		same, length = boot.find_id(translated[0], translated[1])
# 		identity = float(same) / float(length)


# 		if identity < 0.5:
# 			if to_write[k[0]][-8:] == 'AAAAAAAA':
# 				to_write[k[1]] = mark.complementary(to_write[k[1]])
# 			else: 
# 				to_write[k[0]] = mark.complementary(to_write[k[0]])

# 			infodict = cut.framecut_only_straight(to_write, '/mnt/gamma/user/sofya/scripts/final/dnds/2_coding_parts_unaligned.fasta'\
# 				, front, end, start)

# 			if infodict[k[0]] == 'not_coding' or infodict[k[1]] == 'not_coding':
# 				outline += '\tnc\tnc\tnc\tnc\tnc\tnc\tnc'
# 				continue 

# 			line = ['perl', '/mnt/gamma/user/sofya/tools/translatorx.pl', '-i', \
# 			'/mnt/gamma/user/sofya/scripts/final/dnds/2_coding_parts_unaligned.fasta',\
# 			 '-o', '/mnt/gamma/user/sofya/scripts/final/dnds/3_aligned_seqs', \
# 			 '-p', 'M', '-t', 'F', '-w', '1', '-c', '10']	
# 			process = subprocess.Popen(line)
# 			process.communicate()
# 			# pause = int(input())

# 			translated = reader.readfasta_tolist('/mnt/gamma/user/sofya/scripts/final/dnds/3_aligned_seqs.nt_ali.fasta')
# 			sd, S, nd, N = ds.getdnds(translated[0], translated[1])

# 			same, length = boot.find_id(translated[0], translated[1])
# 			identity = float(same) / float(length)
# 			if identity < 0.5 :
# 				with open('stranges.txt', 'a') as strange:
# 					strange.write(filename + '\t' + k[0] + '\t' + k[1] + '\n')
# 					strange.write(filename + '\t' + k[0] + '\t' + k[1] + '\n')

# 		try: 
# 			pnps = float(nd)/ float(N) / float(sd) * float(S)
# 		except ZeroDivisionError:
# 			pnps = 'NULL'
# #------------------------------------------ align end ----------------------------------------
# 		outline1 += '\t' + '\t'.join(list(map(str, [len(to_write[k[0]]), sd, S, nd, N, pnps, same, length, identity])))
# 		outline2 += '\t' + '\t'.join(list(map(str, [len(to_write[k[0]]), sd, S, nd, N, pnps, same, length, identity])))

# 		print( infodict)
# 		print(translated)
# 		print(filename)
# 		# pause = int(input())
		
# 		for index in range(len(infodict[k[0]])):
# 			somestop = infodict[k[0]][index]
# 			if index == 0:
# 				before = 0
# 			else:
# 				before = infodict[k[0]][index - 1]
# 			piece_len = somestop - before

# 			al_before, al_somestop = get_inalignmment(translated[0], before, somestop)

# 			sd, S, nd, N = ds.getdnds(translated[0][al_before:al_somestop], translated[1][al_before:al_somestop])
# 			same, length = boot.find_id(translated[0][al_before:al_somestop], translated[1][al_before:al_somestop])			
# 			try: 
# 				pnps = float(nd)/ float(N) / float(sd) * float(S)
# 			except ZeroDivisionError:
# 				pnps = 'NULL'

# 			outline1 += '\t' + '\t'.join(list(map(str, [piece_len, sd, S, nd, N, pnps, same, length, identity])))

# 		for index in range(len(infodict[k[1]])):
# 			somestop = infodict[k[1]][index]
# 			if index == 0:
# 				before = 0
# 			else:
# 				before = infodict[k[1]][index - 1]
# 			piece_len = somestop - before

# 			al_before, al_somestop = get_inalignmment(translated[1], before, somestop)

# 			sd, S, nd, N = ds.getdnds(translated[0][al_before:al_somestop], translated[1][al_before:al_somestop])
# 			try: 
# 				pnps = float(nd)/ float(N) / float(sd) * float(S)
# 			except ZeroDivisionError:
# 				pnps = 'NULL'
# 			same, length = boot.find_id(translated[0][al_before:al_somestop], translated[1][al_before:al_somestop])
# 			outline2 += '\t' + '\t'.join(list(map(str, [piece_len, sd, S, nd, N, pnps, same, length, identity])))
		


# info.close()