def get_inalignmment(ald_seq, posstart, posend):
	# print(posstart, posend)
	alstart = 0 
	alend = 0 
	filled = 0 
	now = -1 
	
	while ((filled <= posstart) and (now < len(ald_seq))) :
		now += 1
		if now == len(ald_seq):
			break
		try: 
			if ald_seq[now] != '-':
				filled += 1
		except IndexError:
			print(now)
			print(len(allseq))
			print(posend)
			print(posstart)

			pause = int(input())
	
	alstart = now
	# print(alstart)

	while ((filled <= posend) and  (now < len(ald_seq))) :
		now += 1
		if now == len(ald_seq):
			break
		try: 
			if ald_seq[now] != '-':
				filled += 1
		except IndexError:
			print('ERROR IN GET_INALIGNMENT')
			print(now)
			print(len(allseq))
			print(posend)
			print(posstart)
			pause = int(input())

	alend = now
	# print(alstart, alend)
	return alstart, alend





def clean(seq):
	new_seq = str()
	for elem in seq:
		if elem != '-':
			new_seq += elem
	return new_seq

def position_align_to_norm(alt_seq, position):
	filled = 0 
	now = -1 
	while  (now < position):
		now += 1
		if now > len(alt_seq):
			print("ERROR:: alignment_library.py :: function position_align_to_norm() :: position  > alignment.length")
			pause = int(input())
		if alt_seq[now] != '-':
			filled += 1
	return filled


def position_norm_to_align(alt_seq, position):
	filled = -1 
	now = -1 
	while  (filled < position):
		now += 1
		if now > len(alt_seq):
			print("ERROR:: alignment_library.py :: function position_norm_to_align() :: position  > alignment.length")
			pause = int(input())
		if alt_seq[now] != '-':
			filled += 1
	return now

print(norm_alignment('-AAs--Assa', 2))