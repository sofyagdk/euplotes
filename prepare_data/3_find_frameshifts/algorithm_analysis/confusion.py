from optparse import OptionParser
import sys
sys.path.append("../../../lib/")
import reader

# def FalseNegative(borders, seq, window):
#     FN = 0
#     if borders[0] > window:
#         if sum(seq[:borders[0] - window]) > 0: #seq[border[0] - window:border[0]]
#             FN += 1
#     if len(seq) - borders[1] > window:
#         if sum(seq[borders[1] + window :]) > 0: #seq[border[1] :border[1] + window]
#             FN += 1
#     return FN



def FalseNegative5(borders, seq, window):
    FN = 0
    if borders[0] > window:
        if sum(seq[:borders[0] - window]) > 0: #seq[border[0] - window:border[0]]
            FN += 1
    return FN


def FalseNegative3(borders, seq, window):
    FN = 0
    if len(seq) - borders[1] > window:
        if sum(seq[borders[1] + window :]) > 0: #seq[border[1] :border[1] + window]
            FN += 1
    return FN

def get_bp_seq(idr):
    with open('idr_to_filename.txt', 'r') as fin:
        for line in fin:
            copy = line.split()
            if copy[0] == idr:
                filename = copy[1]
                break
    path = '../../data/nuc_alignments/' + filename
    bp_dict = reader.readfasta_todictionary(path)
    bp_dict = reader.filter_paralogues(bp_dict)

    return bp_dict[idr]



def scan(seq, coord):
    laststop = coord
    while ((laststop > 3) and ((seq[laststop:laststop + 3] in ['TAA', 'TAG']) == False)):
        laststop -= 3
    return laststop


def FalsePositive(coding, main, seq, seq_bp):
    unknown = 0 
    control_Reads = 0 
    control_length = 0 

    FP5, FP_length5, FP3, FP_length3 = 0, 0, 0, 0
    TP5, TP3 = 0, 0
    
    for i in range(1, len(coding), 2):
        if coding[i] == main[1]:
            control_Reads = sum(seq[coding[i-1]:coding[i]])
            control_length = coding[i] - coding[i-1]
        
        elif control_length == 0: # meaning we havent yet been in main frame
            # m = min(coding[i], main[0])
            m = scan(seq_bp, coding[i+1])
            if (((m - coding[i-1]) > 15 ) and ( sum(seq[coding[i-1]:m - 15]) == 0)):
                FP5 += 1
                FP_length5 = coding[i] - coding[i-1]
            elif ((coding[i] - coding[i-1]) <= 15):
                # print('IM HERE')
                unknown += 1
                # pause = int(input())
            elif (((m - coding[i-1]) > 15 ) and ( sum(seq[coding[i-1]:m - 15]) != 0)):
                TP5 += 1


        elif control_length != 0:
            if (sum(seq[coding[i-1] + 40:coding[i]]) == 0):
                FP3 += 1
                FP_length3 = coding[i] - coding[i-1]
            else:
                TP3 += 1

    return control_Reads//30, control_length, FP_length5,  FP_length3, FP5, FP3, TP5, TP3, unknown



def get_coding(idr, codingsfilename):
    with open(codingsfilename, 'r') as fin:
        for line in fin:
            copy = line.split()
            #hom_cras_c9739_g1_i1.fasta minu_c24042_g1_i1 87,402,404,1307 272,1307 87,402,404,1307 272,1307
            if copy[1] == idr:
                coding = copy[2].split(',')
                main = copy[3].split(',')
                return list(map(int, coding)), list(map(int,  main))




parser = OptionParser()
parser.add_option("-f", "--filename", help="")
parser.add_option("-c", "--code", help="")
opt, args = parser.parse_args()


codingsfilename  =  opt.filename
code = opt.code
# codingsfilename = '../10_parametertests/models/utlength35_20lt_and_IDfit.txt'

FN_3 = [0]*100
FN_5 = [0]*100
falseposfile = open('./falses/false_positive_4_'+ codingsfilename.split('/')[-1]  , 'w')
falseposfile.write('FP5now FP_length5 FP3now FP_length3 control_Reads control_length\n')
TP5 = 0 
FP5 = 0 
TP3 = 0 
FP3 = 0 
allinall = 0 

UN5 = 0
doubles = 0

fs = 0 
unfound = list()
with open('./riboseq/cras_arch.txt', 'r') as fin:
    for line in fin:
        if line[0] == '>':
            allinall += 1
            info = line[1:-1]
            info = info.split()
            coding = info[3].split(',')  
            coding = list(map(int, coding))
            try:
                coding, main = get_coding(info[0], codingsfilename)
                seq_bp = get_bp_seq(info[0])
                seq_bp = reader.clean(seq_bp)
            except:
                print(line)
                print('ho')
                print(unfound)
                unfound.append(info[0])
                continue
        else:
            seq = line[:-1].split()
            seq = list(map(int, seq))
            for window in range(100):
                FN_3[window] += FalseNegative3([coding[0], coding[-1]], seq, window)
                FN_5[window] += FalseNegative5([coding[0], coding[-1]], seq, window)

            if len(coding) < 3:
                continue

            fs += len(coding)//2 - 1
            control_Reads, control_length, FP_length5, FP_length3, FP5now, FP3now, TP5now, TP3now, UN5now = FalsePositive(coding, main, seq, seq_bp)
            if TP5now > 1:
                doubles += 1
            FP5 += FP5now
            TP5 += TP5now
            FP3 += FP3now
            TP3 += TP3now
            UN5 += UN5now
            if (FP5now + FP3now) == 0:
                continue
            else:
                falseposfile.write(' '.join(list(map(str, [FP5now, FP_length5, FP3now, FP_length3, control_Reads, control_length]))) + '\n')
falseposfile.close()

print(fs)
with open('./falses/false_negative_4_'+ codingsfilename.split('/')[-1] , 'w') as fout:
    for window in range(100):
        fout.write(str(window) + ' ' + str(FN_3[window])+ ' ' + str(FN_5[window]) + '\n')


with open('unfound.txt', 'w') as fout:
    fout.write(' '.join(unfound) + '\n')


print('doubles: ', doubles)



FN5 =  FN_5[14]

TN5 = allinall - FP5 - FN5 - TP5
A5  = FP5/float(TN5 + FP5)
B5 = FN5/float(TP5 + FN5)
print('must be same', FN5, TP5,  TP5 + FN5)



FN3 = FN_3[39]
TN3 = allinall - FP3 - FN3 - TP3
A3  = FP3/float(TN3 + FP3)
B3 = FN3/float(TP3 + FN3)
print('must be same', FN3, TP3, TP3 + FN3)

print('False positive:', FP3, FP5)
print('True positive:', TP3, TP5)
print(allinall)

with open('models_parameter_vfin.txt', 'a') as fout:
    fout.write(' '.join(list(map(str, [codingsfilename, TN5, FP5, FN5, TP5, A5, B5, UN5, TN3, FP3, FN3, TP3, A3, B3, code ]))) + '\n')