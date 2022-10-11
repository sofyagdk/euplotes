from ete3 import Tree


class String_Tree(Tree):

	def set_sequence(self, sequence):
		self.sequence = sequence

	def set_shifts(self, shiftlist):

		self.shift_events = list(set(shiftlist))
		self.shift_events.sort()

	def add_shift(self, position):
		self.shift_events.append(position)
		self.shift_events.sort()

	def tree_distance(self):
		d = 0
		for child in self._children:
			d += child.tree_distance()
			d += child.dist
		return d

	def count_pos_mutations(self, pos, mutto):

		root_nuc = self.sequence[pos]
		if mutto == root_nuc:
			return 0
		mutnum = 0

		for child in self.children:
			if child.sequence[pos] == mutto:
				mutnum += 1
			elif child.sequence[pos] == root_nuc:
				mutnum += child.count_pos_mutations(pos, mutto)

		return mutnum



	# def get_mutations(self):

	# 	answer = list()


	# 	for node in self.traverse("postorder"): 
	# 		if node.is_leaf():
	# 			continue
	# 		if node.is_root():
	# 			continue

	# 		else:
	# 			for codon in seq:
	# 				if codon == lys:

	# 					lysine_type = 

	# 					for child:
	# 						if child_codon == lys:
	# 							branch_length = 






    # preorder: 1)Visit the root, 2) Traverse the left subtree , 3) Traverse the right subtree.
    # postorder: 1) Traverse the left subtree , 2) Traverse the right subtree, 3) Visit the root
    # levelorder (default): every node on a level before is visited going to a lower level


	def solve_shifts(self):
		Sure_dict = dict()
		Question_dict = dict()
		self.Event_dict = dict() # position:[ancestor_node, child_node, type = {'gain', 'loss'}]


		for node in self.traverse("postorder"): 
			### sure - more than a half ancestors have the event
			### question - A half ancsestors have the event

			if node.is_leaf():
				Sure_dict[node.name] = list(node.shift_events)
				Question_dict[node.name] = list()
				continue
		

			ch1 = node.children[0].name
			ch2 = node.children[1].name
			positions = list(set(Sure_dict[ch1]).union(Sure_dict[ch2], Question_dict[ch1], Question_dict[ch2]))

			Sure_dict[node.name] = list()
			Question_dict[node.name] = list()

			for pos in positions:
				if (pos in Sure_dict[ch1]) and (pos in Sure_dict[ch2]):
					Sure_dict[node.name].append(pos)
				elif ((pos in Sure_dict[ch1]) or (pos in Sure_dict[ch2])) and ((pos in Question_dict[ch1]) or (pos in Question_dict[ch2])):
					Sure_dict[node.name].append(pos)
				elif (pos in Question_dict[ch1]) and (pos in Question_dict[ch2]):
					Question_dict[node.name].append(pos)
				elif (pos in Sure_dict[ch1]) or (pos in Sure_dict[ch2]):
					Question_dict[node.name].append(pos)


			Sure_dict[node.name] = list(set(Sure_dict[node.name]))
			Question_dict[node.name] = list(set(Question_dict[node.name]))


		for node in self.traverse("preorder"):
			if node.is_root():
				node.shift_events = list(set(Sure_dict[node.name]))
				node.shift_events.sort()
				continue

			ancestor_shifts = node.up.shift_events
			ancestor_name = node.up.name

			
			for anc_pos in ancestor_shifts:
				if (anc_pos in Sure_dict[node.name]) or (anc_pos in Question_dict[node.name]):
					node.shift_events.append(anc_pos)
				else:
					if (pos in self.Event_dict) == False:
						self.Event_dict[pos] = list()
					self.Event_dict[pos].append(ancestor_name)
					self.Event_dict[pos].append(node.name)
					self.Event_dict[pos].append('loss')

			
			for node_pos in Sure_dict[node.name]:

				if node_pos in ancestor_shifts:
					continue
				else:
					node.shift_events.append(node_pos)
					if (node_pos not in self.Event_dict):
						self.Event_dict[node_pos] = list()
					self.Event_dict[node_pos].append(ancestor_name)
					self.Event_dict[node_pos].append(node.name)
					self.Event_dict[node_pos].append('gain')

			node.shift_events = list(set(node.shift_events))
			
			# print node.name
