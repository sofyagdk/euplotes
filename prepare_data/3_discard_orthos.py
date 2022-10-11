




import os 
import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
import reader

### very high ID
# to_disc = ['eury_eury_c28021_g1_i1.fasta',
#  'eury_eury_c27830_g1_i1.fasta',
#  'minu_c2325_g2_i1.fasta',
#  'minu_c35148_g1_i1.fasta',
#  'eury_eury_c20509_g1_i1.fasta',
#  'minu_c37094_g1_i1.fasta',
#  'eury_eury_c22606_g1_i1.fasta',
#  'minu_c34808_g1_i1.fasta',
#  'eury_eury_c23894_g1_i1.fasta',
#  'minu_c38898_g1_i1.fasta',
#  'eury_eury_c8798_g1_i1.fasta',
#  'minu_c23449_g1_i1.fasta',
#  'eury_eury_c29393_g1_i1.fasta',
#  'eury_eury_c25709_g1_i1.fasta',
#  'eury_eury_c28542_g1_i1.fasta',
#  'eury_eury_c3407_g1_i1.fasta',
#  'minu_c31461_g1_i1.fasta',
#  'eury_eury_c28165_g1_i1.fasta',
#  'minu_c6620_g1_i1.fasta',
#  'minu_c4765_g1_i1.fasta',
#  'eury_eury_c25431_g1_i1.fasta',
#  'foca_c30695_g1_i1.fasta',
#  'cras_c5896_g1_i1++.fasta',
#  'minu_c33317_g1_i1.fasta',
#  'minu_c28790_g1_i1.fasta',
#  'minu_c25493_g1_i1.fasta',
#  'cras_c32031_g1_i1++.fasta',
#  'eury_eury_c28438_g1_i1.fasta',
#  'eury_eury_c3436_g1_i1.fasta',
#  'eury_eury_c12172_g1_i1.fasta',
#  'cras_c3363_g1_i1++.fasta']


# for elem in to_disc:
# 	os.system('cp ./nuc_alignments/' + elem + ' ./ortho_discarded/' + elem)
# 	os.system('rm ./nuc_alignments/' + elem)





low_id_and_pair = ['minu_c30453_g1_i1.fasta', 
'harp_c35937_g1_i1.fasta', 
'eury_eury_c27051_g1_i1.fasta', 
'minu_c27279_g1_i1.fasta', 
'eury_eury_c24871_g1_i1.fasta', 
'cras_c5325_g1_i1++.fasta', 
'cras_c2265_g1_i1++.fasta', 
'eury_eury_c19165_g1_i1.fasta', 
'eury_eury_c28832_g1_i1.fasta', 
'eury_eury_c2845_g1_i1.fasta', 
'eury_eury_c25499_g1_i1.fasta', 
'eury_eury_c28919_g1_i1.fasta', 
'eury_eury_c14622_g1_i1.fasta', 
'harp_c37569_g1_i1.fasta', 
'cras_c9847_g1_i1++.fasta', 
'eury_eury_c9631_g1_i2.fasta', 
'minu_c37853_g1_i1.fasta', 
'harp_c53601_g1_i1.fasta', 
'cras_c8677_g1_i1++.fasta', 
'harp_c16376_g1_i1.fasta', 
'minu_c34574_g1_i1.fasta', 
'eury_eury_c23535_g1_i1.fasta', 
'minu_c9322_g1_i1.fasta', 
'cras_c5207_g1_i2.fasta', 
'eury_eury_c22850_g1_i1.fasta', 
'minu_c30056_g1_i1.fasta', 
'minu_c9887_g2_i1.fasta', 
'eury_eury_c22500_g1_i1.fasta', 
'eury_eury_c23664_g1_i1.fasta', 
'eury_eury_c11278_g1_i1.fasta', 
'minu_c39356_g1_i1.fasta', 
'cras_c4054_g1_i1++.fasta', 
'minu_c28307_g1_i1.fasta', 
'eury_eury_c23914_g1_i1.fasta', 
'minu_c26518_g1_i1.fasta', 
'eury_eury_c11238_g1_i1.fasta', 
'minu_c27322_g1_i1.fasta', 
'cras_c10919_g1_i1++.fasta', 
'eury_eury_c19475_g1_i1.fasta', 
'minu_c26393_g1_i1.fasta', 
'minu_c23743_g1_i1.fasta', 
'minu_c5631_g1_i1.fasta', 
'eury_eury_c9332_g1_i1.fasta', 
'eury_eury_c29579_g1_i1.fasta', 
'foca_c30629_g1_i1.fasta', 
'minu_c32726_g1_i1.fasta', 
'eury_eury_c22397_g1_i1.fasta', 
'minu_c10790_g2_i1.fasta', 
'eury_eury_c23175_g1_i1.fasta', 
'minu_c21770_g1_i1.fasta', 
'eury_eury_c25652_g1_i1.fasta', 
'foca_c2659_g1_i1.fasta', 
'eury_eury_c5391_g1_i1.fasta', 
'eury_eury_c1821_g1_i1.fasta', 
'eury_eury_c4590_g1_i1.fasta', 
'minu_c6774_g1_i1.fasta', 
'cras_c5315_g1_i1++.fasta', 
'harp_c61494_g1_i1.fasta', 
'petz_comp27857_c0_seq1.fasta', 
'eury_eury_c5349_g1_i1.fasta']






# for elem in low_id_and_pair:
# 	os.system('cp ./nuc_alignments/' + elem + ' ./ortho_discarded/' + elem)
# 	os.system('rm ./nuc_alignments/' + elem)



# to_turn = {'cras_c5007_g1_i1++.fasta' : 'eury_eury_c24506_g1_i1',
# 'cras_c2913_g1_i1++.fasta' : 'harp_c36544_g1_i1',
# 'harp_c37141_g1_i1.fasta' : 'octo_c59502_g1_i1',
# 'cras_c13999_g1_i1++.fasta' : 'minu_c33616_g1_i1',
# 'eury_eury_c12004_g1_i1.fasta' : 'rari_c12411_g1_i1',
# 'cras_c3983_g1_i1++.fasta' : 'cras_c3983_g1_i1++',
# 'eury_eury_c7859_g1_i1.fasta' : 'rari_c4970_g1_i1',
# 'cras_c24471_g1_i1++.fasta' : 'raik_c20881_g4_i1',
# 'cras_c26090_g1_i1++.fasta' : 'harp_c863_g1_i1',
# 'harp_c38223_g2_i1.fasta' : 'octo_c19351_g1_i1',
# 'cras_c7213_g1_i1++.fasta' : 'minu_c17356_g1_i1',
# 'cras_c5286_g1_i1++.fasta' : 'minu_c38662_g1_i1',
# 'cras_c16907_g1_i1++.fasta' : 'cras_c16907_g1_i1++',
# 'cras_c29940_g1_i1++.fasta' : 'minu_c10983_g1_i1',
# 'cras_c3177_g1_i1++.fasta' : 'minu_c27710_g1_i1',
# 'cras_c6551_g1_i1++.fasta' : 'harp_c37942_g1_i1',
# 'eury_eury_c18763_g1_i1.fasta' : 'eury_eury_c18763_g1_i1',
# 'cras_c1667_g1_i1++.fasta' : 'foca_c30119_g1_i1',
# 'cras_c6179_g1_i1++.fasta' : 'eury_eury_c9241_g1_i1',
# 'cras_c31199_g1_i1++.fasta' : 'cras_c31199_g1_i1++',
# 'cras_c19612_g1_i1.fasta' : 'rari_c3777_g1_i1',
# 'cras_c26380_g1_i1++.fasta' : 'rari_c5500_g1_i1',
# 'eury_eury_c417_g1_i1.fasta' : 'rari_c2239_g1_i1',
# 'eury_eury_c28398_g1_i1.fasta' : 'rari_c7180_g1_i1',
# 'cras_c17093_g1_i1++.fasta' : 'cras_c17093_g1_i1++',
# 'cras_c5595_g1_i1++.fasta' : 'rari_c6359_g1_i1',
# 'cras_c29525_g1_i1+cras_comp7339_c0_seq1+.fasta' : 'minu_c37546_g1_i1',
# 'cras_c8091_g1_i1++.fasta' : 'cras_c8091_g1_i1++',
# 'eury_eury_c15759_g1_i1.fasta' : 'petz_comp26527_c0_seq1',
# 'eury_eury_c28419_g1_i1.fasta' : 'harp_c35858_g1_i1',
# 'cras_c2497_g1_i1++.fasta' : 'harp_c37864_g2_i1'}


# compliment = {'A':'T', 'T':'A', 'G':'C', 'C':'G', '-':'-'}
# def reverse_compliment(seq):
# 	revseq = list()
# 	for letter in seq:
# 		revseq.append(compliment[letter])
# 	return ''.join(revseq[::-1])


# def makefasta_fromdict(bp_dict, fastaname):
# 	with open(fastaname, 'w') as fasta:
# 		for idr in bp_dict:
# 			fasta.write('>' + idr + '\n')
# 			fasta.write(bp_dict[idr] + '\n')

# def align(fasta, alignmentfilename):
# 	cmd = 'muscle -in ' + fasta + ' -out ' + alignmentfilename #also folder name
# 	os.system(cmd)




# for filename in to_turn:
# 	path = './nuc_alignments/'+filename
# 	bp_dict = reader.readfasta_todictionary(path)
# 	idr_to_turn = to_turn[filename]
# 	bp_dict[idr_to_turn] = reverse_compliment(bp_dict[idr_to_turn])

# 	makefasta_fromdict(bp_dict, path)
	
# 	align(path, path)


# manualy_discarded = ['eury_eury_c6111_g1_i1.fasta',
# 'cras_c11411_g1_i1++.fasta', 
# 'cras_c26070_g1_i1++.fasta',
# 'eury_eury_c21575_g1_i1.fasta']

# manualy_truned = ['eury_eury_c6111_g1_i1.fasta', 
# 'eury_eury_c3474_g1_i1.fasta', 
# 'cras_c4123_g1_i1++.fasta', 
# 'cras_c2207_g1_i1++.fasta',
# 'cras_c6295_g1_i1++.fasta', 
# 'cras_c6353_g1_i1++.fasta', 
# 'cras_c8085_g1_i1++.fasta']  




too_high_id = ['eury_eury_c20351_g1_i1.fasta',
 'eury_eury_c23982_g1_i1.fasta',
 'cras_c31199_g1_i1++.fasta']


to_low_id = ['eury_eury_c19056_g1_i1.fasta',
 'eury_eury_c21411_g1_i1.fasta',
 'eury_eury_c23417_g1_i1.fasta',
 'eury_eury_c24564_g1_i1.fasta',
 'eury_eury_c27252_g1_i1.fasta',
 'eury_eury_c27478_g1_i1.fasta',
 'eury_eury_c29036_g1_i1.fasta',
 'eury_eury_c4590_g1_i2.fasta',
 'foca_c10387_g2_i1.fasta',
 'raik_c41277_g1_i1.fasta',
 'cras_c22558_g1_i1++.fasta',
 'minu_c384_g1_i1.fasta',
 'cras_comp14802_c0_seq1.fasta',
 'minu_c24735_g1_i1.fasta',
 'minu_c23467_g1_i1.fasta',
 'minu_c35102_g1_i1.fasta']


for elem in to_low_id:
	os.system('cp ./nuc_alignments/' + elem + ' ./ortho_discarded/' + elem)
	os.system('rm ./nuc_alignments/' + elem)


for elem in too_high_id:
	os.system('cp ./nuc_alignments/' + elem + ' ./ortho_discarded/' + elem)
	os.system('rm ./nuc_alignments/' + elem)


