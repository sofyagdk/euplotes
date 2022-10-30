# sftp://sofya@mg.uncb.iitp.ru/mnt/gamma/user/sofya/data/Riboseq/E.crassus_Ribo.fq

# def eprint(*args, **kwargs):
#     print(*args, file=sys.stderr, **kwargs)




pool = [0]*9620943
out = dict()
multicopy = dict()

with open('./ribototrans.tsv', 'r') as fin:
	for line in fin:
		copy = line.split()
		idr = copy[1]
		
		# if (('cras' in copy[1] ) == False):
			# continue
			# print(line[:-1]) # '''python kart.py > otheranimal.tsv'''
		
		read = int(copy[0])
		
		if (pool[read] != 0):
			pool[read] = pool[read] + ',' + idr
			multicopy[read] = pool[read] 
			# print(last[:-1]) 
			# print(line[:-1]) # '''python kart.py > manytimes.tsv '''

		pool[read] = idr

		if (('cras' in copy[1] ) == False):
			continue
		if idr not in out:
			out[idr] = list(['', ''])
		start = int(copy[8])
		end = int(copy[9])
		if end > start: # normal
			out[idr][0] += str(start) + ',' + str(end) + ','
		else:
			out[idr][1] += str(start) + ',' + str(end) + ','

for idr in out:
	outline = '\t'.join([idr, str(len(out[idr][0])//2), str(len(out[idr][1])//2), out[idr][0], out[idr][1]])
	print(outline)

with open('mapcrassus_light.txt', 'w') as fout:
	for idr in out:
		outline = '\t'.join([idr, str(len(out[idr][0])//2), str(len(out[idr][1])//2)])
		fout.write(outline + '\n')

# with open('multicopy.txt', 'w') as fout:
# 	for read in multicopy:
# 		fout.write(str(read) + '\t' + multicopy[read] + '\n')
#
# with open('multicopy_nocrass.txt', 'w') as fout:
# 	for read in multicopy:
# 		if 'cras' not in multicopy[read]:
# 			fout.write(str(read) + '\t' + multicopy[read] + '\n')


# 7866175	raik_c21408_g2_i1	100.000	28	0	0	1	28	999	1026	6.20e-09	52.8
# 		start = int(copy[8])
# 		end = int(copy[9])
# 		if end < start:
# 			print(line[:-1])
# 			revids.add( copy[1] )

# 		last = line
# # print('finite first part')
# print(len(revids))


# bits = 0 
# with open('./blastwork/ribobits.fasta', 'r') as fin:
# 	for line in fin:
# 		if line[0]  ==  '>':
# 			bits += 1

# eprint(bits)
# eprint(sum(pool))
# eprint(sum(pool)/float(bits))



# 9620943
# 585108
# 0.0608160759294












# directory = '../nuc_alignments/'
# files = os.listdir(directory)
# cras = 0
# for filename in files:

# 	path = directory + filename
# 	bp_dict = reader.readfasta_todictionary(path)
# 	bp_dict = reader.filter_paralogues(bp_dict)
# 	idrs = bp_dict.keys()

# 	for idr in bp_dict:
# 		if 'cras' in idr:
# 			cras += 1
# # 		cut, has_shift = getcut(bp_dict, idr, filename)
# print(cras)


# "makeblastdb -in ./blastwork/wholedata.fasta -parse_seqids -out ./blastwork/transcripts_db -dbtype 'nucl' -hash_index"

# blastn -db ./blastwork/transcripts_db -query ./blastwork/example_ofribobits.fasta -outfmt 7 -out ./blastwork/example.tsv -perc_identity 80
# -parse_seqids

# makeblastdb -in test.fsa -parse_seqids -blastdb_version 5 -taxid_map test_map.txt -title "Cookbook demo" -dbtype prot