####################################
# sofya 20.01.2022
# sofya.gaydukova@gmail.com

####################################


import os 

directory = './contamination_fastas/'
filenames = os.listdir(directory)

for filename in filenames:
	dbfilename = filename
	dbfilename.replace('.fna_nt', '')
	dbfilename.replace('.fna', '')
	dbfilename += '_db'
	
	cmd = 'makeblastdb -in ' + directory + filename + ' -parse_seqids -blastdb_version 5 -out ./contatmination_dbs/' + dbfilename + ' -dbtype nucl'
	os.system(cmd)

	with open('db_list.txt', 'a') as fout:
		fout.write(dbfilename + '\n')

		