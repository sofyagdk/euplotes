####################################
# sofya 20.01.2022
# sofya.gaydukova@gmail.com

####################################

import os

#####################################
## INPUT
#####################################

jobname = 'filter_'


db_names = list()
with open('db_list.txt', 'r') as fin:
	for line in fin:
		db_names.append( line.strip()[:-3])


result_dir= '/results/'


tsvdir_temp = '/tsvs_temp/'
tsv_files = os.listdir(tsvdir_temp)

#####################################
## FILLING OUT idr_match_dict
#####################################


# hits = open(result_dir + jobname + 'blast_hits_big_table.tsv', 'w')

idr_match_dict = dict()

with open('../../data/concatenated_transcriptomes.fasta', 'r') as fin:
	for line in fin:
		if (line[0] == '>'):
			idr = line.strip().split()[0][1:]
			idr_match_dict[idr] = list()




for filename in tsv_files:
	path = tsvdir_temp + filename
	currdb = filename[7:-4]

	with open(path, 'r') as fin:
		for line in fin:
			if line[0] == '#':
				continue
			else:
				idr = line.split()[0]
				# if (idr not in idr_match_dict):
				# 	idr_match_dict[idr] = list()
				idr_match_dict[idr].append(currdb)
				# print("--------")
				# print(idr_match_dict[idr]) 
				# pause = int(input())
				# hits.write(line)
				# rari_c23877_g1_i1	NC_050098.1	79.325	237	45	4	237	471	216307489	216307255	9.09e-39	163

# hits.close()


#####################################
## CLEANING - SURE WANNA DO IT?
#####################################

#####################################
# for filename in tsv_files:
# 	path = tsvdir + filename
# 	cmd = 'rm ' + path
# 	os.system(cmd)
#####################################


#####################################
## WRITING RESULTS
#####################################

with open(result_dir + jobname + 'blast_table.txt' ,'w') as fout:
	fout.write('idr ')
	fout.write(' '.join(db_names))
	fout.write('\n')

	for idr in idr_match_dict:
		out = [idr]
		for db in db_names:
			if db in idr_match_dict[idr]:
				out.append('1')
			else:
				out.append('0')
		fout.write(' '.join(out))
		fout.write('\n')

