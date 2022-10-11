# deletions     AAATAA   4   100
# insertions  AAA-A[AG]  35  1076
# insertions  AAG-A[AG]  15  3540
import pandas as pd

def get_special_context(table, context_list_before, central_letter,  context_list_after):
    returner = list()
    i = 0 
    for row in table.itertuples():
        if row[2][5 - len(context_list_before[0]):5] in context_list_before and row[2][6:6+len(context_list_after[0])] in context_list_after and (row[2][5] in central_letter or row[3][5] in central_letter):
            returner.append(True)
        else:
            returner.append(False)
                
        i+=1
        
    return table.iloc[returner]

def get_special_letter(table, central_letter):
    returner = list()
    i = 0 
    for row in table.itertuples():
        if (row[2][5] in central_letter or row[3][5] in central_letter):
            returner.append(True)
        else:
            returner.append(False)
                
        i+=1
        
    return table.iloc[returner]


res = open('shi_treehairbrush.txt', 'w')
neu = open('neu_treehairbrush.txt', 'w')



events_all = pd.read_csv('./events_treehairbrush_golden.txt', header = 0, delim_whitespace=True)

returner = list()
i = 0 
for row in events_all.itertuples():
	if '-' in row[2][:5]  or '-' in row[2][6:]:
		returner.append(False)
		continue
	else:
		returner.append(True)
	i+=1
events = events_all.iloc[returner]

# res.write('deletions\tAAATAA\t' + str(len(events[events['event_type'] == 'deletion'][events['place'] == 'shi'])) + '\t' + str(94) + '\n')
res.write('deletions\tAAATAA\t' + str(len(events[events['event_type'] == 'deletion'][events['place'] == 'shi'])) + '\t' + str(104) + '\n')



returner = list()
i = 0 
for row in events.itertuples():
    if row[2][4]  == 'G':
        returner.append(True)
        continue
    else:
        returner.append(False)
    i+=1
events_g = events.iloc[returner]

res.write('insertions\tAAG-A[AG]\t' + str(len(events_g[events_g['event_type'] == 'insertion'][events_g['place'] == 'shi'])) + '\t'  ) #+ str(94) + '\n')

with open('contexts_forced_treehairbrush_gold.txt', 'r') as context:
	for line in context:
		if 'AAGA[GAR][CAGTRYSWKMBDHVN]' in line:
			copy = line.split()
			res.write(copy[1] + '\n')



returner = list()
i = 0 
for row in events.itertuples():
    if row[2][4]  == 'A':
        returner.append(True)
        continue
    else:
        returner.append(False)
    i+=1
events_a = events.iloc[returner]
res.write('insertions\tAAA-A[AG]\t' + str(len(events_g[events_g['event_type'] == 'insertion'][events_g['place'] == 'shi'])) + '\t'  ) #+ str(94) + '\n')

with open('contexts_forced_treehairbrush_gold.txt', 'r') as context:
	for line in context:
		if 'AAAA[GAR][CAGTRYSWKMBDHVN]' in line:
			copy = line.split()
			res.write(copy[1] + '\n')




neu1 = get_special_context(events, ["AA"], ["T"],  ['A'])
neu1 = neu1[neu1['place'] == 'ins']
ad = len(neu1[neu1['event_type'] == 'deletion'])
ai = len(neu1[neu1['event_type'] == 'insertion'])
neu1 = get_special_context(events, ["AG"], ["T"],  ['A'])
neu1 = neu1[neu1['place'] == 'ins']
gd = len(neu1[neu1['event_type'] == 'deletion'])
gi = len(neu1[neu1['event_type'] == 'insertion'])

de = 0

with open('contexts_forced_treehairbrush_gold.txt', 'r') as context:
	for line in context:
		if '[CAGTRYSWKMBDHVN]AGTA[CAGTRYSWKMBDHVN]' in line:
			copy = line.split()
			de += int(copy[1])
		if '[CAGTRYSWKMBDHVN]AATA[CAGTRYSWKMBDHVN]' in line:
			copy = line.split()
			de += int(copy[1])

neu.write('deletions\tAATA\t' + str(ad + gd) + '\t' + str(de) + '\n')


with open('contexts_forced_treehairbrush_gold.txt', 'r') as context:
	for line in context:
		if '[CAGTRYSWKMBDHVN]AAA[CAGTRYSWKMBDHVN][CAGTRYSWKMBDHVN]' in line:
			copy = line.split()
			neu.write('insertions\tAATA\t' + str(ai) + '\t'+ copy[1] + '\n')
with open('contexts_forced_treehairbrush_gold.txt', 'r') as context:
	for line in context:
		if '[CAGTRYSWKMBDHVN]AGA[CAGTRYSWKMBDHVN][CAGTRYSWKMBDHVN]' in line:
			copy = line.split()
			neu.write('insertions\tAGTA\t' + str(gi) + '\t' + copy[1] + '\n')



res.close()
neu.close()


# res = open('shi_tree2.txt', 'w')
# neu = open('neu_tree2.txt', 'w')



# events_all = pd.read_csv('./events_tree2_final.txt', header = 0, delim_whitespace=True)

# returner = list()
# i = 0 
# for row in events_all.itertuples():
# 	if '-' in row[2][:5]  or '-' in row[2][6:]:
# 		returner.append(False)
# 		continue
# 	else:
# 		returner.append(True)
# 	i+=1
# events = events_all.iloc[returner]

# # res.write('deletions\tAAATAA\t' + str(len(events[events['event_type'] == 'deletion'][events['place'] == 'shi'])) + '\t' + str(94) + '\n')
# res.write('deletions\tAAATAA\t' + str(len(events[events['event_type'] == 'deletion'][events['place'] == 'shi'])) + '\t' + str(72) + '\n')



# returner = list()
# i = 0 
# for row in events.itertuples():
#     if row[2][4]  == 'G':
#         returner.append(True)
#         continue
#     else:
#         returner.append(False)
#     i+=1
# events_g = events.iloc[returner]

# res.write('insertions\tAAG-A[AG]\t' + str(len(events_g[events_g['event_type'] == 'insertion'][events_g['place'] == 'shi'])) + '\t'  ) #+ str(94) + '\n')

# with open('contexts_forced_tree2_final.txt', 'r') as context:
# 	for line in context:
# 		if 'AAGA[GAR][CAGTRYSWKMBDHVN]' in line:
# 			copy = line.split()
# 			res.write(copy[1] + '\n')



# returner = list()
# i = 0 
# for row in events.itertuples():
#     if row[2][4]  == 'A':
#         returner.append(True)
#         continue
#     else:
#         returner.append(False)
#     i+=1
# events_a = events.iloc[returner]
# res.write('insertions\tAAA-A[AG]\t' + str(len(events_g[events_g['event_type'] == 'insertion'][events_g['place'] == 'shi'])) + '\t'  ) #+ str(94) + '\n')

# with open('contexts_forced_tree2_final.txt', 'r') as context:
# 	for line in context:
# 		if 'AAAA[GAR][CAGTRYSWKMBDHVN]' in line:
# 			copy = line.split()
# 			res.write(copy[1] + '\n')




# neu1 = get_special_context(events, ["AA"], ["T"],  ['A'])
# neu1 = neu1[neu1['place'] == 'ins']
# ad = len(neu1[neu1['event_type'] == 'deletion'])
# ai = len(neu1[neu1['event_type'] == 'insertion'])
# neu1 = get_special_context(events, ["AG"], ["T"],  ['A'])
# neu1 = neu1[neu1['place'] == 'ins']
# gd = len(neu1[neu1['event_type'] == 'deletion'])
# gi = len(neu1[neu1['event_type'] == 'insertion'])

# de = 0

# with open('contexts_forced_tree2_final.txt', 'r') as context:
# 	for line in context:
# 		if '[CAGTRYSWKMBDHVN]AGTA[CAGTRYSWKMBDHVN]' in line:
# 			copy = line.split()
# 			de += int(copy[1])
# 		if '[CAGTRYSWKMBDHVN]AATA[CAGTRYSWKMBDHVN]' in line:
# 			copy = line.split()
# 			de += int(copy[1])

# neu.write('deletions\tAATA\t' + str(ad + gd) + '\t' + str(de) + '\n')


# with open('contexts_forced_tree2_final.txt', 'r') as context:
# 	for line in context:
# 		if '[CAGTRYSWKMBDHVN]AAA[CAGTRYSWKMBDHVN][CAGTRYSWKMBDHVN]' in line:
# 			copy = line.split()
# 			neu.write('insertions\tAATA\t' + str(ai) + '\t'+ copy[1] + '\n')
# with open('contexts_forced_tree2_final.txt', 'r') as context:
# 	for line in context:
# 		if '[CAGTRYSWKMBDHVN]AGA[CAGTRYSWKMBDHVN][CAGTRYSWKMBDHVN]' in line:
# 			copy = line.split()
# 			neu.write('insertions\tAGTA\t' + str(gi) + '\t' + copy[1] + '\n')



# res.close()
# neu.close()