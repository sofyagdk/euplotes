####################################
# sofya 20.01.2022
# sofya.gaydukova@gmail.com

# program that writes cmd into bash file 
# blast_against_comtamination.sh
# if the resulting blast sh file doesnt work -- maybe the full path needs to be
# specified (line 34)
# It is also posssible to change jobname
####################################

query_sequences_file = '../../data/concatenated_transcriptomes.fasta'
jobname = 'filter_'


dbs = list()
db_dir = './contatmination_dbs/'

with open('db_list.txt', 'r') as fin:
	for line in fin:
		dbs.append(db_dir + line.strip())

fout = open('blast_against_comtamination.sh', 'w') 
fout.write('#!/bin/bash\n#$ -cwd\n')


for db_name in dbs:
	cmd = 'blastn -db ' + db_name
	cmd += ' -query ' + query_sequences_file + ' -outfmt 7'
	db_clearname = db_name.split('/')[-1]
	cmd += ' -out ./tsvs_temp/' + jobname + db_clearname[:-3] + '.tsv' ### here you may need a full path to file
	cmd += ' -evalue 1e-25 -perc_identity 70'
	fout.write(cmd + '\n')


