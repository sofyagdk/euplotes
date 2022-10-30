import sys
sys.path.append("../../../../lib/")
import reader
import os


def SpeciesIdrInDictionary(keys, sp):
    for index, idr in enumerate(keys):
        if sp == idr[:4]:
            return idr
    return None

directory = '../../../data/nuc_alignments/'
fout = open('crassus_transcripts.fasta', 'w')
for filename in os.listdir(directory):
    path = directory +filename
    bp_dict = reader.readfasta_todictionary(path)
    idr = SpeciesIdrInDictionary(bp_dict.keys(), 'cras')
    if idr:
        fout.write('>' + idr[:50] + '\n')
        fout.write(bp_dict[idr].replace('-', '')+ '\n')

fout.close()