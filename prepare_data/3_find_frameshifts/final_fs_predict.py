import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader
import ds
import random 
from optparse import OptionParser
import numpy as np
import os
import math
import pandas as pd
import scipy.stats

file_idr = pd.read_csv('../file_to_idr.txt', header = None, delim_whitespace=True)


logsfilename ='../logs_graph.txt'


logs = pd.read_csv(logsfilename, header = None, delim_whitespace=True)

logs.rename(columns = {0:'idr', 
1:'start', #in alignment
2:'end', # in alignment
3:'length', # in bp
4:'tomain', # main/3/5
5:'is3', #true/false
6:'biggest_prot_ID', # only float value
7:'best_prot_ID', # id,same,ungapped - and bigget ungapped
8:'all_id_data', # id,same,ungapped x as many times as calculated
9:'fits_ID', # true/false
10:'bestfit_prot_ID', 
11:'best_pnps', # for main
12:'best_counting_list', # for main
13:'fit_best_pnps', # a float from 0 to 1  - fiiting_value
14:'pnps_fit_list', # all calculated fiiting_value s
15:'all_pnps' # dnds,nd,N,sd,S,fitting_value x as many times as calculated
}, inplace = True)
# 16:'best_paml_w',
# 17:'best_paml_data',
# 18:'fit_best_paml_w',
# 19:'paml_w_fit_list',
# 20:'all_paml_data'}, inplace = True)


fitlook = list()
whigh = list()
for index, row in logs.iterrows():
    if  type(row['all_pnps']) is float:
    	fitlook.append('None')
    	whigh.append('None')
        continue

    e = list(map(float, row['all_pnps'].split(',')))
    b_N = 0
    v = 1
    w = 1
    for i in range(0, len(e), 6):
        dnds, nd, N, sd, S, st = e[i:i + 6]
        if b_N < N:
            b_N = N
            v = st
            w = dnds
    fitlook.append(v)
    whigh.append(w)
logs['fithigh'] = fitlook   
logs['whigh'] = whigh  


# logs = logs[logs['is3'] == False]


def alignment_norm(alt_seq, position):
	filled = -1 
	now = -1 
	while  (now < position):
		now += 1
		if now >= len(alt_seq):
			break
		if alt_seq[now] != '-':
			filled += 1
	return filled



def  alignment_norm_coding(alt_seq, coding):
	norm_coding = list()
	for pos in coding:
		newpos = alignment_norm(alt_seq, pos)
		norm_coding.append(newpos)
	return norm_coding
# def is_frame_coding(row):
# 	if row['fits_ID'] == 'True':
# 		return True
# 	return False


# ut5 = 100
# ut3 = 100

# lt5 = 20
# lt3 = 36

# fout = open('./models/lengths_' + str(lt5) + '-' + str(ut5) + '_' + str(lt3) + '-' + str(ut3) + '.txt', 'w')
idmin = 0.7
x = 0 

def is_frame_coding(row, ut5, ut3, lt5, lt3, s, e, alt_seq):



	if ((row['start'] < s) and (row['end'] > e)):
		return False


	if row['best_prot_ID'] == 'None':
		return False
	else:
		identity, same, ungapped = row['best_prot_ID'].split(',')
		identity = float(identity)

	dndsfit = 0.85
	# if ((row['fit_best_pnps'] == 'None') or (float(row['fit_best_pnps']) > dndsfit)):
	# 	return False

	
	if ((row['start'] < s)  and (row['end'] > s)):

		n1, n2 = alignment_norm_coding(alt_seq, [row['start'], s])#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		# n1, n2 = alignment_norm_coding(alt_seq, [row['start'], row['end']])
		# if (n2 - n1)>= 3*ut5:
			# return True


		if (n2 - n1) < 3*lt5 :
			return False
		elif (n2 - n1) >= 3*ut5:
			if (identity > 0.64):
				# if ((row['whigh'] == 'None') or (float(row['whigh']) > 0.5)):
				# 	return False
				return True 
		else:
			if row['fits_ID'] in ['True', True]:
				if ((row['fithigh'] == 'None') or (float(row['fithigh']) > dndsfit)):
					return False
				return True 
			


		return False



	if (row['end'] > e):
		# if row['tomain'] !='3':

		

		n1, n2 = alignment_norm_coding(alt_seq, [e, row['end']])
 

		if (n2 - n1) < 3*lt3 :
			return False
		elif (n2 - n1) >= 3*ut3:
			if (identity > 0.64):
				# if ((row['whigh'] == 'None') or (float(row['whigh']) > 0.5)):
				# 	return False
				return True 
		else:
			if row['fits_ID'] in ['True', True]:
				if ((row['fithigh'] == 'None') or (float(row['fithigh']) > dndsfit)):
					return False
				return True 
			


		return False



	
	return False

# def is_frame_coding(row, ut, s, e, alt_seq):
# 	# if row['fits_ID'] == 'True':
# 	# 	return True
# 	if ((row['length'] >= ut*30) and (row['start'] < s )):
# 		return True
# 	if (row['end'] > e ):
		
# 		n1, n2 = alignment_norm_coding(alt_seq, [e, row['end']])
# 		if (n2 - n1)>= 3*ut:
# 			return True 
# 	return False




def get_mainali(idr):
	temp = logs[logs['idr'] == idr]
	mainali = list([int(temp[temp['tomain'] == 'main']['start']), int(temp[temp['tomain'] == 'main']['end'])])
	return mainali


def norm_alignment(alt_seq, position):
	filled = -1 
	now = -1 
	while  (filled < position):
		now += 1
		if alt_seq[now] != '-':
			filled += 1
	return now


def norm_alignment_coding(alt_seq, coding):
	al_coding = list()
	for pos in coding:
		newpos = norm_alignment(alt_seq, pos)
		al_coding.append(newpos )
	return al_coding

def trans(coding):
	for i in range(2, len(coding), 2):
		nextone = coding[i - 1] + 1
		if nextone + 2 < coding[i]:
			print('NONONO')
			break 
		coding[i] = nextone 
		while (coding[i+1] - coding[i])%3 != 0:
			coding[i] += 1



	return coding



file_idr = pd.read_csv('../file_to_idr.txt', header = None, delim_whitespace=True)
# file_idr

def get_alignment(myidr):
	filename= file_idr[file_idr[0] == myidr][1]
	path =  '../data/nuc_alignments/' + filename
	bp_dict = reader.readfasta_todictionary(path)
	bp_dict = reader.filter_paralogues(bp_dict)

	# for idr in bp_dict:
	# 	if idr == myidr:
	return bp_dict[idr]







def get_framesali(idr, ut5, ut3, lt5, lt3, alt_seq):
	# global logs
	temp = logs[logs['idr'] == idr]
	# print(temp)

	if len(temp[temp['tomain'] == 'main']['start']) > 1:
		print(len(temp[temp['tomain'] == 'main']['start']) , idr)

	s = int(list(temp[temp['tomain'] == 'main']['start'])[0])
	e = int(list(temp[temp['tomain'] == 'main']['end'])[0])
	main = alignment_norm_coding(alt_seq, [s, e])

	# framesali.append([s, e])
	# print(frames)
	
#     print(temp)
	temp = temp[temp['length'] >= min(lt3, lt5)*3]
	temp = temp[temp['tomain'] != 'main']


	fittingframes = list()
	fittingframes.append([s,e])


	for i in range(7):
		s =  fittingframes[0][0]
		e = fittingframes[-1][1]
	

		for index, row in temp.iterrows():
			# if is_frame_coding(row, ut, s, e, alt_seq):
			if (is_frame_coding(row, ut5, ut3, lt5, lt3, s, e, alt_seq) == True):
				fr_s = row['start']
				fr_e = row['end']
				if ((fr_s < fittingframes[-1][1]) and (fr_e > fittingframes[-1][1])): # from 3'
					fittingframes.append([fr_s, fr_e])
				elif ((fr_s < fittingframes[0][0]) and (fr_e > fittingframes[0][0])): #from 5'
					bef = [[fr_s, fr_e]]
					bef.extend(fittingframes)
					fittingframes = list(bef)
			s =  fittingframes[0][0]
			e = fittingframes[-1][1]
	
	# return frames




# def turn_aliframes_to_coding(framesali, mainali, alt_seq, ut5, ut3, lt5, lt3):
	# frames = list()

	# s = int(temp[temp['tomain'] == 'main']['start'])
	# e = int(temp[temp['tomain'] == 'main']['end'])

	frames_norm = list()

	for frame in fittingframes:
		newfr = alignment_norm_coding(alt_seq, frame)
		frames_norm.append(newfr)

	# fittingframes = list()
	# fittingframes.append(main)

	# for i in range(7):

	# 	for fr in frames:
	# 		if fr == main:
	# 			continue
	# 		if fr[0] < fittingframes[-1][1] and fr[1] > fittingframes[-1][1]: # from 3'
	# 			if fr[1] - fittingframes[-1][1] > lt3*3:
	# 				if is_frame_coding(row, ut5, ut3, lt5, lt3, s, e, alt_seq)
	# 				fittingframes.append(fr)

	# 		if fr[0] < fittingframes[0][0] and fr[1] > fittingframes[0][0]: #from 5'
	# 			if (fr[1] - fr[0] > lt5*3) and (fittingframes[0][1] - fr[1] > lt5*3):
	# 				bef = list([fr])
	# 				bef.extend(fittingframes)
	# 				fittingframes = list(bef)

	coding = list()
	for frame in frames_norm:
		coding.extend(frame)



	
	coding = trans(coding)

	al_coding = norm_alignment_coding(alt_seq, coding)
	al_main = norm_alignment_coding(alt_seq, main)

	return coding, main, al_coding, al_main






directory = '/mnt/gamma/user/sofya/scripts/final/0pipeline/nuc_alignments/'
files  = os.listdir(directory)



# def is_frame_coding(row, ut5, ut3, lt5, lt3, s, e, alt_seq): # redefinition here!


#  	global idmaxn
# 	if row['best_prot_ID'] == 'None':
# 		return False
# 	else:
# 		identity, same, ungapped = row['best_prot_ID'].split(',')
# 		identity = float(identity)

# 	if ((row['start'] < s) and (row['end'] > e)):
# 		return False

# 	elif (row['start'] < s):
# 		n1, n2 = alignment_norm_coding(alt_seq, [row['start'], s])
# 		if (n2 - n1)>= 3*ut5:
# 			return True
# 		if (n2 - n1)< 3*lt5:
# 			return False

# 		if identity > 0.85:
# 			return True

# 	elif (row['end'] > e):
# 		n1, n2 = alignment_norm_coding(alt_seq, [e, row['end']])
# 		if ((n2 - n1) >= 3*ut3):
# 			return True 
# 		else:
# 			return False
	
# 	return False



ut3 = 60
lt3  = 25

ut5 = 60
lt5 = 15

fout = open('./models/final_7.txt', 'w')
justfilename = 'final_7.txt'


for filename in files:
	path = directory + filename
	bp_dict = reader.readfasta_todictionary(path)
	bp_dict = reader.filter_paralogues(bp_dict)

	for idr in bp_dict:
		coding, main, coding_ali, main_ali = get_framesali(idr, ut5, ut3, lt5, lt3, bp_dict[idr])
		fout.write(filename + ' ' + idr + ' ' + ','.join(list(map(str, coding)))+ ' ' + ','.join(list(map(str, [main[0],main[1]]))) + ' ' + ','.join(list(map(str, coding_ali)))+' ' + ','.join(list(map(str, main_ali)))+ '\n')

fout.close()

with open('../9_riboseq/conf_pipe_7.sh', 'a') as fout:
	fout.write('python confusion.py -f ../10_parametertests/models/' + justfilename + ' -c fin  \n'  )

