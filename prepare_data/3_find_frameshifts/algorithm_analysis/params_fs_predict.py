import sys

sys.path.append("../../../lib/")
import reader
import os
import pandas as pd

logsfilename = '../logs_graph.txt'

logs = pd.read_csv(logsfilename, header=None, delim_whitespace=True)

logs.rename(columns={0: 'idr',
                     1: 'start',  # in alignment
                     2: 'end',  # in alignment
                     3: 'length',  # in bp
                     4: 'tomain',  # main/3/5
                     5: 'is3',  # true/false
                     6: 'biggest_prot_ID',  # only float value
                     7: 'biggest_ungapped_prot_ID',  # id,same,ungapped - and bigget ungapped
                     8: 'all_id_data',  # id,same,ungapped x as many times as calculated
                     9: 'fits_ID',  # true/false
                     10: 'bestfit_prot_ID',
                     11: 'best_pnps',  # for main
                     12: 'best_counting_list',  # for main
                     13: 'fit_best_pnps',  # a float from 0 to 1  - fitting_value
                     14: 'pnps_fit_list',  # all calculated fiiting_value s
                     15: 'all_pnps'  # dnds,nd,N,sd,S,fitting_value x as many times as calculated
                     }, inplace=True)


# logs = logs[logs['is3'] == False]

# 16:'best_paml_w',
# 17:'best_paml_data',
# 18:'fit_best_paml_w',
# 19:'paml_w_fit_list',
# 20:'all_paml_data'}, inplace = True)

def alignment_norm(alt_seq, position):

    filled = -1
    now = -1
    while (now < position):
        now += 1
        if now >= len(alt_seq):
            return filled
        if alt_seq[now] != '-':
            filled += 1
    return filled


def alignment_norm_coding(alt_seq, coding):
    norm_coding = list()
    for pos in coding:
        newpos = alignment_norm(alt_seq, pos)
        norm_coding.append(newpos)
    return norm_coding


def get_mainali(idr):
    temp = logs[logs['idr'] == idr]
    mainali = list([int(temp[temp['tomain'] == 'main']['start']), int(temp[temp['tomain'] == 'main']['end'])])
    return mainali


def norm_alignment(alt_seq, position):
    filled = -1
    now = -1
    while (filled < position):
        now += 1
        if alt_seq[now] != '-':
            filled += 1
    return now


def norm_alignment_coding(alt_seq, coding):
    al_coding = list()
    for pos in coding:
        newpos = norm_alignment(alt_seq, pos)
        al_coding.append(newpos)
    return al_coding


def trans(coding):
    for i in range(2, len(coding), 2):
        nextone = coding[i - 1] + 1
        if nextone + 2 < coding[i]:
            print('NONONO')
            break
        coding[i] = nextone
        while (coding[i + 1] - coding[i]) % 3 != 0:
            coding[i] += 1

    return coding


def get_fitted_pnps(row):
    separated = row['all_pnps'].split(',')
    fit_value = float(row['fit_best_pnps'])

    for index, elem in enumerate(separated[::6]):

        pnps_list = separated[index * 6: (index + 1) * 6]
        if float(pnps_list[-1]) == fit_value:
            return float(pnps_list[0])

    return None


# fout = open('./models/lengths_' + str(hard_length_5) + '-' + str(soft_length_5) + '_' + str(hard_length_3) + '-' + str(soft_length_3) + '.txt', 'w')
# idmin = 0.7
# x = 0


def check_short_frame(row,
                      fits_identity=True,
                      min_identity=0.4,
                      dnds_uniform=0.85,
                      dnds_hard_val=0):

    if fits_identity:
        if row['fits_ID'] == 'False':
            return False

    if min_identity:
        if row['biggest_prot_ID'] == 'None':
            return False
        if float(row['biggest_prot_ID']) < min_identity:
            return False

    if dnds_uniform:
        if row['fit_best_pnps'] == 'None' or float(row['fit_best_pnps']) < dnds_uniform:
            return False

    if dnds_hard_val:
        dnds = get_fitted_pnps(row)
        if dnds < dnds_hard_val:
            return False

    return True


def check_long_frame(row, min_identity_for_long):
    if row['biggest_ungapped_prot_ID'] == 'None':
        return False
    else:
        identity, same, ungapped = row['biggest_ungapped_prot_ID'].split(',')
        identity = float(identity)

    if min_identity_for_long:
        if identity > min_identity_for_long:
            return True

    if float(row['biggest_prot_ID']) > identity:
        print('Found frame where ungapped ID is too low but the biggest ID is high')
        print(row['all_id_data'])

    return False


def is_frame_coding(row, s, e, alt_seq,
                    soft_length_3=60, soft_length_5=60,
                    hard_length_3=15, hard_length_5=25,
                    fits_identity=True, min_identity=0.4,
                    min_identity_for_long=0.64,
                    dnds_uniform=0.85,
                    dnds_hard_val=0):
    if (row['start'] < s) and (row['end'] > e):
        return False

    if ((row['start'] < s) and (row['end'] > s)):
        n1, n2 = alignment_norm_coding(alt_seq, [row['start'], s])

        if (n2 - n1) < 3 * hard_length_5:
            return False

        if (n2 - n1) >= 3 * soft_length_5:
            return check_long_frame(row, min_identity_for_long)

        return check_short_frame(row, fits_identity, min_identity, dnds_uniform, dnds_hard_val)

    if (row['end'] > e):

        n1, n2 = alignment_norm_coding(alt_seq, [e, row['end']])

        if (n2 - n1) < 3 * hard_length_3:
            return False

        if (n2 - n1) >= 3 * soft_length_3:
            return check_long_frame(row, min_identity_for_long)

        return check_short_frame(row, fits_identity, min_identity, dnds_uniform, dnds_hard_val)

    return False


idr_to_filename = pd.read_csv('idr_to_filename.txt', header=None, delim_whitespace=True)


def get_alignment(myidr):
    filename = idr_to_filename[idr_to_filename[0] == myidr][1]
    path = '../../data/nuc_alignments/' + filename
    bp_dict = reader.readfasta_todictionary(path)
    bp_dict = reader.filter_paralogues(bp_dict)

    return bp_dict[idr]


def get_coding_frames(idr, alt_seq,
                      soft_length_3=60, soft_length_5=60,
                      hard_length_3=15, hard_length_5=25,
                      fits_identity=True, min_identity=0.4,
                      min_identity_for_long=0.64,
                      dnds_uniform=0.85,
                      dnds_hard_val=0):
    global logs
    temp = logs[logs['idr'] == idr]

    print(idr)

    s = int(list(temp[temp['tomain'] == 'main']['start'])[0])
    e = int(list(temp[temp['tomain'] == 'main']['end'])[0])


    main = alignment_norm_coding(alt_seq, [s, e])

    temp = temp[temp['length'] >= min(hard_length_3, hard_length_5) * 3]
    temp = temp[temp['tomain'] != 'main']

    fittingframes = list()
    fittingframes.append([s, e])

    for i in range(15):
        s = fittingframes[0][0]
        e = fittingframes[-1][1]

        for index, row in temp.iterrows():
            if (is_frame_coding(row, s, e, alt_seq,
                                soft_length_3, soft_length_5,
                                hard_length_3, hard_length_5,
                                fits_identity, min_identity, min_identity_for_long,
                                dnds_uniform, dnds_hard_val) == True):

                fr_s = row['start']
                fr_e = row['end']
                if ((fr_s < fittingframes[-1][1]) and (fr_e > fittingframes[-1][1])):  # from 3'
                    fittingframes.append([fr_s, fr_e])
                elif ((fr_s < fittingframes[0][0]) and (fr_e > fittingframes[0][0])):  # from 5'
                    bef = [[fr_s, fr_e]]
                    bef.extend(fittingframes)
                    fittingframes = list(bef)
            s = fittingframes[0][0]
            e = fittingframes[-1][1]

    frames_norm = list()

    for frame in fittingframes:
        newfr = alignment_norm_coding(alt_seq, frame)
        frames_norm.append(newfr)

    coding = list()
    for frame in frames_norm:
        coding.extend(frame)

    coding = trans(coding)

    al_coding = norm_alignment_coding(alt_seq, coding)
    al_main = norm_alignment_coding(alt_seq, main)

    return coding, main, al_coding, al_main


directory = '../../data/nuc_alignments/'
files = os.listdir(directory)

########################################################
# Block 1. Hard thresholds.
########################################################

#################### Upper threshold for 3' ####################################


# soft_length_3 = 1000
# hard_length_3 = 1000
# soft_length_5 = 1000
# hard_length_5 = 1000
#
# # for hard_length_3 in range(10, 155, 5):
# for hard_length_5 in [5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24] + list(range(25, 155, 5)):
#     # soft_length_5 = hard_length_5
#     # soft_length_3 = hard_length_3
#
#     justfilename = 'lengths_' + str(hard_length_5) + '-' + str(soft_length_5) + \
# 				   '_' + str(hard_length_3) + '-' + str(soft_length_3) + '.txt'
#
#     fout = open('./models/' + justfilename, 'w')
#
#     for filename in files:
#         path = directory + filename
#         bp_dict = reader.readfasta_todictionary(path)
#         bp_dict = reader.filter_paralogues(bp_dict)
#
#         for idr in bp_dict:
#             coding, main, coding_ali, main_ali = get_coding_frames(idr, bp_dict[idr], soft_length_3, soft_length_5,hard_length_3, hard_length_5,
#                                                                    fits_identity = False, min_identity = 0, min_identity_for_long = 0, dnds_uniform = 0, dnds_hard_val = 0)
#
#
#             fout.write(filename + ' ' + idr + ' ' + ','.join(list(map(str, coding))) + ' ' + ','.join(
#                 list(map(str, [main[0], main[1]]))) + ' ' + ','.join(list(map(str, coding_ali))) + ' ' + ','.join(
#                 list(map(str, main_ali))) + '\n')
#
#     fout.close()
#
#     with open('confusion_pipe.sh', 'a') as fout:
#         fout.write('python confusion.py -f ./models/' + justfilename + ' -c hard_length_3  \n')

########################################################

########################################################
# Block 2. Soft thresholds.
########################################################


soft_length_3 = 1000
hard_length_3 = 25
soft_length_5 = 1000
hard_length_5 = 15

# for hard_length_3 in range(10, 155, 5):
for soft_length_3 in range(hard_length_3 + 5, 155, 5):
# for soft_length_5 in range(hard_length_5, 155, 5):

    justfilename = 'lengths_' + str(hard_length_5) + '-' + str(soft_length_5) + \
				   '_' + str(hard_length_3) + '-' + str(soft_length_3) + '.txt'

    fout = open('./models/' + justfilename, 'w')

    for filename in files:
        path = directory + filename
        bp_dict = reader.readfasta_todictionary(path)
        bp_dict = reader.filter_paralogues(bp_dict)

        for idr in bp_dict:
            coding, main, coding_ali, main_ali = get_coding_frames(idr, bp_dict[idr], soft_length_3, soft_length_5, hard_length_3, hard_length_5,
                                                                   fits_identity = False, min_identity = 0, min_identity_for_long = 0, dnds_uniform = 0, dnds_hard_val = 0)


            fout.write(filename + ' ' + idr + ' ' + ','.join(list(map(str, coding))) + ' ' + ','.join(
                list(map(str, [main[0], main[1]]))) + ' ' + ','.join(list(map(str, coding_ali))) + ' ' + ','.join(
                list(map(str, main_ali))) + '\n')

    fout.close()

    with open('confusion_pipe.sh', 'a') as fout:
        fout.write('python confusion.py -f ./models/' + justfilename + ' -c soft_length_3  \n')








########################################################
# Block 3. Minimal identity
########################################################

# soft_length_5 = 1000
# hard_length_5 = 15
# soft_length_3 = 1000
# hard_length_3 = 25
#
# # for hard_length_3 in range(10, 155, 5):
# # for soft_length_3 in range(hard_length_3, 155, 5):
# for min_id in [0, 0.1, 0.2, 0.3,  0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
#
#     justfilename = 'lengths_15-1000-25-1000_' +'min_id_' +str(min_id)+ '.txt'
#
#     fout = open('./models/' + justfilename, 'w')
#
#     for filename in files:
#         path = directory + filename
#         bp_dict = reader.readfasta_todictionary(path)
#         bp_dict = reader.filter_paralogues(bp_dict)
#
#         for idr in bp_dict:
#             coding, main, coding_ali, main_ali = get_coding_frames(idr, bp_dict[idr], soft_length_3, soft_length_5, hard_length_3, hard_length_5,
#                                                                    fits_identity=False, min_identity = min_id, min_identity_for_long=min_id, dnds_uniform=0, dnds_hard_val=0)
#
#
#             fout.write(filename + ' ' + idr + ' ' + ','.join(list(map(str, coding))) + ' ' + ','.join(
#                 list(map(str, [main[0], main[1]]))) + ' ' + ','.join(list(map(str, coding_ali))) + ' ' + ','.join(
#                 list(map(str, main_ali))) + '\n')
#
#     fout.close()
#
#     with open('confusion_pipe.sh', 'a') as fout:
#         fout.write('python confusion.py -f ./models/' + justfilename + ' -c min_id  \n')
#











# for soft_length_5 in range(10, 155, 5):
# 	hard_length_5  = soft_length_5 - 1

# 	soft_length_3 = 1000
# 	hard_length_3 = 1000


# 	fout = open('./models/upperlengths_for5ADD_' + str(hard_length_5) + '-' + str(soft_length_5) + '_' + str(hard_length_3) + '-' + str(soft_length_3) + '.txt', 'w')
# 	justfilename = 'upperlengths_for5ADD_' + str(hard_length_5) + '-' + str(soft_length_5) + '_' + str(hard_length_3) + '-' + str(soft_length_3) + '.txt'


# 	for filename in files:
# 		path = directory + filename
# 		bp_dict = reader.readfasta_todictionary(path)
# 		bp_dict = reader.filter_paralogues(bp_dict)

# 		for idr in bp_dict:
# 			coding, main, coding_ali, main_ali = get_coding_frames(idr, soft_length_5, soft_length_3, hard_length_5, hard_length_3, bp_dict[idr])
# 			fout.write(filename + ' ' + idr + ' ' + ','.join(list(map(str, coding)))+ ' ' + ','.join(list(map(str, [main[0],main[1]]))) + ' ' + ','.join(list(map(str, coding_ali)))+' ' + ','.join(list(map(str, main_ali)))+ '\n')

# 	fout.close()

# 	with open('../9_riboseq/conf_pipe_5.sh', 'a') as fout:
# 		fout.write('python confusion.py -f ../10_parametertests/models/' + justfilename + ' -c soft_length_5ADD  \n'  )


####################### Threshold for 5' frame is being calculated from the end of frame #################################

# def is_frame_coding(row, soft_length_5, soft_length_3, hard_length_5, hard_length_3, s, e, alt_seq): # redefinition here!

# 	if ((row['start'] < s) and (row['end'] > e)):
# 		return False

# 	elif (row['start'] < s):
# 		n1, n2 = alignment_norm_coding(alt_seq, [row['start'], row['end']]) # HERE is the difference with the previous
# 		if (n2 - n1)>= 3*soft_length_5:
# 			return True
# 		else:
# 			return False

# 	elif (row['end'] > e):
# 		n1, n2 = alignment_norm_coding(alt_seq, [e, row['end']])
# 		if ((n2 - n1) >= 3*soft_length_3):
# 			return True
# 		else:
# 			return False

# 	return False


# for soft_length_5 in range(10, 155, 5):
# 	hard_length_5  = soft_length_5 - 1

# 	soft_length_3 = 1000
# 	hard_length_3 = 1000


# 	fout = open('./models/upperlengths_for5_' + str(hard_length_5) + '-' + str(soft_length_5) + '_' + str(hard_length_3) + '-' + str(soft_length_3) + '.txt', 'w')
# 	justfilename = 'upperlengths_for5_' + str(hard_length_5) + '-' + str(soft_length_5) + '_' + str(hard_length_3) + '-' + str(soft_length_3) + '.txt'


# 	for filename in files:
# 		path = directory + filename
# 		bp_dict = reader.readfasta_todictionary(path)
# 		bp_dict = reader.filter_paralogues(bp_dict)

# 		for idr in bp_dict:
# 			coding, main, coding_ali, main_ali = get_coding_frames(idr, soft_length_5, soft_length_3, hard_length_5, hard_length_3, bp_dict[idr])
# 			fout.write(filename + ' ' + idr + ' ' + ','.join(list(map(str, coding)))+ ' ' + ','.join(list(map(str, [main[0],main[1]]))) + ' ' + ','.join(list(map(str, coding_ali)))+' ' + ','.join(list(map(str, main_ali)))+ '\n')

# 	fout.close()

# 	with open('../9_riboseq/conf_pipe_5.sh', 'a') as fout:
# 		fout.write('python confusion.py -f ../10_parametertests/models/' + justfilename + ' -c soft_length_5  \n'  )


########################################################
# End block 1
########################################################


#######################################################
# BLOCK 2. Maximum identity
#######################################################
# idmaxn = 0.45
#
# def is_frame_coding(row, soft_length_5, soft_length_3, hard_length_5, hard_length_3, s, e, alt_seq): # redefinition here!
# 	global idmaxn
# 	if row['biggest_ungapped_prot_ID'] == 'None':
# 		return False
# 	else:
# 		identity, same, ungapped = row['biggest_ungapped_prot_ID'].split(',')
# 		identity = float(identity)
#
#
# 	if ((row['start'] < s) and (row['end'] > e)):
# 		return False
#
# 	elif (row['start'] < s):
# 		n1, n2 = alignment_norm_coding(alt_seq, [row['start'], s]) # We are using ADD
# 		if (n2 - n1)>= 3*soft_length_5:
# 			return False ############## to estimate freely
# 		elif  (n2 - n1) < 3*hard_length_5:
# 			return False
# 		elif ungapped < hard_length_5:
# 			return False
#
# 	elif (row['end'] > e):
# 		n1, n2 = alignment_norm_coding(alt_seq, [e, row['end']])
# 		if ((n2 - n1) >= 3*soft_length_3):
# 			return False  ############## to estimate freely
# 		elif (n2 - n1) < 3*hard_length_3:
# 			return False
# 		elif ungapped < hard_length_3:
# 			return False
#
#
# 	if identity > idmaxn:
# 		return True
#
#
# 	return False
#
#
#
# soft_length_3 = 60
# soft_length_5 = 60
#
# hard_length_3 = 25
# hard_length_5 = 15
#
#
# idmaxn = 0.45
# while idmaxn < 1:
# 	idmaxn += 0.05
#
# 	fout = open('./models/idmaxn_' + str(idmaxn) +'_dndsfit1_ADD5_' + str(hard_length_5) + '-' + str(soft_length_5) + '_' + str(hard_length_3) + '-' + str(soft_length_3) + '.txt', 'w')
# 	justfilename = 'idmaxn_' + str(idmaxn) +'_dndsfit1_ADD5_' + str(hard_length_5) + '-' + str(soft_length_5) + '_' + str(hard_length_3) + '-' + str(soft_length_3) + '.txt'
#
#
# 	for filename in files:
# 		path = directory + filename
# 		bp_dict = reader.readfasta_todictionary(path)
# 		bp_dict = reader.filter_paralogues(bp_dict)
#
# 		for idr in bp_dict:
# 			coding, main, coding_ali, main_ali = get_coding_frames(idr, soft_length_5, soft_length_3, hard_length_5, hard_length_3, bp_dict[idr])
# 			fout.write(filename + ' ' + idr + ' ' + ','.join(list(map(str, coding)))+ ' ' + ','.join(list(map(str, [main[0],main[1]]))) + ' ' + ','.join(list(map(str, coding_ali)))+' ' + ','.join(list(map(str, main_ali)))+ '\n')
#
# 	fout.close()
#
# 	with open('../9_riboseq/conf_pipe_5.sh', 'a') as fout:
# 		fout.write('python confusion.py -f ../10_parametertests/models/' + justfilename + ' -c idmaxn  \n'  )


#######################################################
# End block 2
#######################################################


#######################################################
# BLOCK 3. Minimum identity
#######################################################
# idmin = 0.4
#
# def is_frame_coding(row, soft_length_5, soft_length_3, hard_length_5, hard_length_3, s, e, alt_seq): # redefinition here!
# 	global idmin
# 	if row['biggest_ungapped_prot_ID'] == 'None':
# 		return False
# 	else:
# 		identity, same, ungapped = row['biggest_ungapped_prot_ID'].split(',')
# 		identity = float(identity)
#
#
# 	if ((row['start'] < s) and (row['end'] > e)):
# 		return False
#
# 	elif (row['start'] < s):
# 		n1, n2 = alignment_norm_coding(alt_seq, [row['start'], s]) # We are using ADD
# 		if (n2 - n1)>= 3*soft_length_5:
# 			return False ############## to estimate freely
# 		elif  (n2 - n1) < 3*hard_length_5:
# 			return False
# 		elif ungapped < hard_length_5:
# 			return False
#
# 	elif (row['end'] > e):
# 		n1, n2 = alignment_norm_coding(alt_seq, [e, row['end']])
# 		if ((n2 - n1) >= 3*soft_length_3):
# 			return False  ############## to estimate freely
# 		elif (n2 - n1) < 3*hard_length_3:
# 			return False
# 		elif ungapped < hard_length_3:
# 			return False
#
#
# 	if ((identity > idmin) and (row['fits_ID'] == True)):
# 		return True
#
#
# 	return False
#
#
#
# soft_length_3 = 60
# soft_length_5 = 60
#
# hard_length_3 = 25
# hard_length_5 = 15
#
#
# idmin = 0.45
# while idmin < 1:
# 	idmin += 0.05
#
# 	fout = open('./models/idmin_' + str(idmin) +'_dndsfit1_ADD5_' + str(hard_length_5) + '-' + str(soft_length_5) + '_' + str(hard_length_3) + '-' + str(soft_length_3) + '.txt', 'w')
# 	justfilename = 'idmin_' + str(idmin) +'_dndsfit1_ADD5_' + str(hard_length_5) + '-' + str(soft_length_5) + '_' + str(hard_length_3) + '-' + str(soft_length_3) + '.txt'
#
#
# 	for filename in files:
# 		path = directory + filename
# 		bp_dict = reader.readfasta_todictionary(path)
# 		bp_dict = reader.filter_paralogues(bp_dict)
#
# 		for idr in bp_dict:
# 			coding, main, coding_ali, main_ali = get_coding_frames(idr, soft_length_5, soft_length_3, hard_length_5, hard_length_3, bp_dict[idr])
# 			fout.write(filename + ' ' + idr + ' ' + ','.join(list(map(str, coding)))+ ' ' + ','.join(list(map(str, [main[0],main[1]]))) + ' ' + ','.join(list(map(str, coding_ali)))+' ' + ','.join(list(map(str, main_ali)))+ '\n')
#
# 	fout.close()
#
# 	with open('../9_riboseq/conf_pipe_5.sh', 'a') as fout:
# 		fout.write('python confusion.py -f ../10_parametertests/models/' + justfilename + ' -c idmin  \n'  )
#

#######################################################
# End block 3
#######################################################


########################################################
# BLOCK 4. DnDs
########################################################
# dndsfit = 0.4
#
# def is_frame_coding(row, soft_length_5, soft_length_3, hard_length_5, hard_length_3, s, e, alt_seq): # redefinition here!
# 	global dndsfit
# 	# if row['biggest_ungapped_prot_ID'] == 'None':
# 	# 	return False
# 	# else:
# 	# 	identity, same, ungapped = row['biggest_ungapped_prot_ID'].split(',')
# 	# 	identity = float(identity)
#
#
# 	if ((row['start'] < s) and (row['end'] > e)):
# 		return False
#
# 	elif (row['start'] < s):
# 		n1, n2 = alignment_norm_coding(alt_seq, [row['start'], s]) # We are using ADD
# 		if (n2 - n1)>= 3*soft_length_5:
# 			return False ############## to estimate freely
# 		elif  (n2 - n1) < 3*hard_length_5:
# 			return False
#
# 	elif (row['end'] > e):
# 		n1, n2 = alignment_norm_coding(alt_seq, [e, row['end']])
# 		if ((n2 - n1) >= 3*soft_length_3):
# 			return False  ############## to estimate freely
# 		elif (n2 - n1) < 3*hard_length_3:
# 			return False
#
# 	if ((row['fit_best_pnps'] != 'None') and (float(row['fit_best_pnps']) < dndsfit)):
# 		return True
#
# 	# if ((identity > idmin) and (row['fits_ID'] == True)):
# 	# 	return True
#
#
# 	return False
#
#
#
# soft_length_3 = 60
# soft_length_5 = 60
#
# hard_length_3 = 25
# hard_length_5 = 15
#
#
# dndsfit = 0.2
# while dndsfit < 1:
# 	dndsfit += 0.05
#
# 	fout = open('./models/dndsfit_' + str(dndsfit) +'_ADD5_' + str(hard_length_5) + '-' + str(soft_length_5) + '_' + str(hard_length_3) + '-' + str(soft_length_3) + '.txt', 'w')
# 	justfilename = 'dndsfit_' + str(dndsfit) +'_ADD5_' + str(hard_length_5) + '-' + str(soft_length_5) + '_' + str(hard_length_3) + '-' + str(soft_length_3) + '.txt'
#
#
# 	for filename in files:
# 		path = directory + filename
# 		bp_dict = reader.readfasta_todictionary(path)
# 		bp_dict = reader.filter_paralogues(bp_dict)
#
# 		for idr in bp_dict:
# 			coding, main, coding_ali, main_ali = get_coding_frames(idr, soft_length_5, soft_length_3, hard_length_5, hard_length_3, bp_dict[idr])
# 			fout.write(filename + ' ' + idr + ' ' + ','.join(list(map(str, coding)))+ ' ' + ','.join(list(map(str, [main[0],main[1]]))) + ' ' + ','.join(list(map(str, coding_ali)))+' ' + ','.join(list(map(str, main_ali)))+ '\n')
#
# 	fout.close()
#
# 	with open('../9_riboseq/conf_pipe_5.sh', 'a') as fout:
# 		fout.write('python confusion.py -f ../10_parametertests/models/' + justfilename + ' -c dndsfit  \n'  )
#

########################################################
# End block 4
########################################################


# soft_length_5 = 65
# soft_length_3 = 60

# # hard_length_5 from 5aa to 65aa
# # hard_length_3 from 14aa to 60aa
# # lets take
# hard_length_5 = 15
# hard_length_3 = 15


# # for soft_length_5 in range(10, 155, 5):
# # 	hard_length_5  = soft_length_5 - 1

# # 	soft_length_3 = 1000
# # 	hard_length_3 = 1000


# fout = open('./models/final.txt' , 'w')
# justfilename  = 'final.txt'


# for filename in files:
# 	path = directory + filename
# 	bp_dict = reader.readfasta_todictionary(path)
# 	bp_dict = reader.filter_paralogues(bp_dict)

# 	for idr in bp_dict:
# 		# print(idr)
# 		# framesali = get_coding_frames(idr, ut, lt, bp_dict[idr])
# 		# framesali = get_coding_frames(idr, soft_length_5, soft_length_3, hard_length_5, hard_length_3, bp_dict[idr])
# 		# print(len(framesali))
# 		# mainali = get_mainali(idr)
# 		# coding, main, coding_ali, main_ali = turn_aliframes_to_coding(framesali, mainali, bp_dict[idr], soft_length_5, soft_length_3, hard_length_5, hard_length_3)
# 		coding, main, coding_ali, main_ali = get_coding_frames(idr, soft_length_5, soft_length_3, hard_length_5, hard_length_3, bp_dict[idr])
# 		# if main_ali != mainali:
# 		# 	print(mainali, main_ali)
# 		# 	print('NONONON')
# 		# 	print(coding, main, coding_ali, main_ali)
# 		# 	pause = int(input())

# 		fout.write(filename + ' ' + idr + ' ' + ','.join(list(map(str, coding)))+ ' ' + ','.join(list(map(str, [main[0],main[1]]))) + ' ' + ','.join(list(map(str, coding_ali)))+' ' + ','.join(list(map(str, main_ali)))+ '\n')

# fout.close()


# """
# unmute!!
# with open('../9_riboseq/conf_pipe_4.sh', 'a') as fout:
# 	fout.write('python confusion.py -f ../10_parametertests/models/' + justfilename + ' -c final  \n'  )
# """
