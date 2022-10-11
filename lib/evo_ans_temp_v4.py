# 20-08-09 sofya added function event_nodes_maker
# 20-04-25 sofya -TEMPORARY ONE
# module for analyzing alignments and finding ancestral state
# 20-07-23 sofya - local module modifyied context decisions

# t = Tree('(minu:0.12310897,cras:0.16946154,(foca:0.26750949,(rari:0.27565106,(eury:0.21540617,((petz:0.28678546,raik:0.26772623)11:0.04916400,(octo:0.34300359,harp:0.21398805)12:0.08092917)13:0.03839098)9:0.04135668)24:0.05101473)30:0.12093409);')
# t = Tree('(petz:0.14339273,(raik:0.26772623,((octo:0.34300359,harp:0.21398805)12:0.08092917,(eury:0.21540617,(rari:0.27565106,(foca:0.26750949,(minu:0.12310897,cras:0.16946154)30:0.12093409)24:0.05101473)9:0.04135668)13:0.03839098)11:0.04916400):0.14339273);')
# print(t.get_tree_root())

# t = Tree('(petz,(raik,((octo,harp),(eury,(rari,(foca,(minu,cras)))))));')

import sys
sys.path.append("/mnt/gamma/user/sofya/scripts/final/chek")
from ete3 import Tree

# coder = dict({['A' , 'G']: 'R', ['C' , 'T'] : 'Y', ['A' , 'T']:'W', ['C', 'G']:'S', ['G' , 'T'] :'K', ['A' , 'C']:'M', ['C' , 'G' , 'T'] : 'B', ['A' , 'G' , 'T'] : 'D', ['A' , 'C' , 'T'] :  'H', ['A' , 'C' , 'G'] : 'V',['A' , 'C' , 'G', 'T'] : 'N'})

def list_to_set(in_list):
	out_set = set()
	for elem in in_list:
		out_set.add(elem)
	return out_set


def to_code(nuc_list):
	if len(nuc_list) == 1:
		return nuc_list[0]

	if len(nuc_list) == 2:
		if ('A' and 'G' in nuc_list):
			return 'R'
		elif ('C' and 'T' in nuc_list):
			return 'Y'
		elif ('G' and 'C' in nuc_list):
			return 'S'
		elif ('A' and 'T' in nuc_list):
			return 'W'
		elif ('G' and 'T' in nuc_list):
			return 'K'
		elif ('A' and 'C' in nuc_list):
			return 'M'

	if len(nuc_list) == 4:
		if 'A' and 'C' and 'G' and 'T' in nuc_list:
			return 'N'

	if len(nuc_list) == 3:
		if ('A' and 'G' and 'T' in nuc_list):
			return 'D'
		elif ('C' and 'G' and 'T' in nuc_list):
			return 'B'
		elif ('A' and 'G' and 'C' in nuc_list):
			return 'V'
		elif ('A' and 'C' and 'T' in nuc_list):
			return 'H'

	return 'X'




	# else:
	# 	nuc_list.sort()
	# 	global coder
	# 	return coder[nuc_list]


def from_code(letter):
	if letter in ['A' , 'C' , 'G', 'T']:
		return letter
	global decoder
	return list(decoder[letter])



class Ancestor_Tree(Tree):

	def set_position(self, pos):
		self.position = pos


	def node_mean(self, mean = None):
		possibles = list()
		for child in self._children:
			if hasattr(child, 'mean'):
				possibles.append(child.mean)
			else:
				child.node_mean()
				possibles.append(child.mean)
		
		if (0 in possibles) and (1 in possibles):
			self.mean = None
		elif (1 in possibles):
			self.mean = 1
		elif (0 in possibles):
			self.mean = 0
		elif self.is_leaf() == True:
			if hasattr(self, 'mean') == False:
				self.mean = mean
		else:
			self.mean = None
	
	def node_mean_downgoing(self, mean = None):
		ancestor_mean = self.mean

		for child in self._children:
			if (child.mean == None) and (child.is_leaf() == False):
				child.mean = ancestor_mean

			child.node_mean_downgoing()

	
	def node_context_upgoing(self, context_length, known_context = None):
		if self.is_leaf() == True:
			if hasattr(self, 'context') == False:
				context = known_context

		context = list()
		possibles = list()

		for child in self._children:
			if hasattr(child, 'context'):
				possibles.append(child.context)
			elif child.is_leaf() == False:
				child.node_context_upgoing(context_length)
				possibles.append(child.context)

		while None in possibles:
			possibles.remove(None)

		# print(possibles)
		# pause = str(input())

		for i in range(context_length):
			letter = set()
			for child_context in possibles:
				letter = letter.union(list_to_set(child_context[i]))

			saver = list_to_set(letter)

			for child_context in possibles:
				letter = letter.intersection(list_to_set(child_context[i]))

			if letter == set():
				letter = list_to_set(saver)

			context.append(list(letter))

		self.context = context


	def node_context_downgoing(self, context_length):
		ancestor_context = self.context

		str_context = str()

		for nuc_list in ancestor_context:
			str_context += to_code(nuc_list)
		
		self.context = str_context

		for child in self._children:
			new_context = list()

			tocheck = child.context
			if tocheck == None:
				continue

			for i in range(context_length):
				inherited = ancestor_context[i]

				inherited_set = list_to_set(inherited)
				tocheck_set = list_to_set(tocheck[i])

				if len(tocheck[i]) == 1:
					new_context.append(list(tocheck[i]))
					continue
				elif len(inherited_set.intersection(tocheck_set)) > 0:
					new_context.append(list(inherited_set.intersection(tocheck_set)))
					continue
				else:
					new_context.append(list(inherited_set.union(tocheck_set)))
					
			child.context = new_context
			child.node_context_downgoing(context_length)

	def tree_distance(self):
		d = 0
		for child in self._children:
			d += child.tree_distance()
			d += child.dist
		return d






	def view_means(self):     ##ASK MISHA
		toprint = self.copy()
		for n in toprint.iter_leaves():
			n.name = n.name +' '+  str(n.mean)
		print(toprint)

	def leaf_context(self, context_string):
		if self.is_leaf() == True:
			self.context = context_string

	def view_contexts(self):
		for n in self.iter_leaves():
			try:
				print(n.context)
			except AttributeError:
				n.context = None
				print(n.context)

	def list_contexts(self):
		context_list = list()
		for n in self.iter_leaves():
			try:
				context_list.append(n.context)
			except AttributeError:
				n.context = None
				context_list.append(n.context)
		return context_list

	def is_onebp_gap(self):
		onebp = False
		for context in self.list_contexts():
			if context == None:
				continue
			if context[len(context)//2] == '-' and context[len(context)//2 - 1] != '-' and context[len(context)//2+1] != '-':
				onebp = True
		return onebp

	
	def make_branch_that_leads_to(self, n):
		after = list()
		before = list()
		for leaf in n.iter_leaves():
			after.append(leaf.name)
		for leaf in self.iter_leaves():
			before.append(leaf.name)	
		before = list(set(before).difference(set(after)))
		branch = [before, after]
		return branch

	def string_branch(self, branch):
		s = str()
		before = branch[0]
		after = branch[1]
		s = ','.join(before)
		s += '_'
		s += ','.join(after)
		return s

	def view_events(self):
		for branch in self.events:
			print(string_branch(self, branch))


	def get_all_events(self):
		meaningful_nodes = list()
		for n in self.iter_leaves():
			if n.mean != None:
				meaningful_nodes.append(n.name)
		meaningful_tree = self.copy()
		meaningful_tree.prune(meaningful_nodes, preserve_branch_length = True)
		meaningful_tree.node_mean()
		meaningful_tree.events_on_branches(meaningful_tree)
		self.events = meaningful_tree.events



	def events_on_branches(self, n):
		if n.is_root() == True:
			self.events = {'deletion': list(), 'insertion': list(), 'unclear': list()}
			for child in self._children:
				self.events_on_branches(child)
		else:
			ancestor_meaning = n._up.mean
			if ancestor_meaning == 1:
				if n.mean == 0:
					self.events['deletion'].append(self.string_branch(self.make_branch_that_leads_to(n)))
				if n.mean == None:
					n.mean = 1

			if ancestor_meaning == 0:
				if n.mean == 1:
					self.events['insertion'].append(self.string_branch(self.make_branch_that_leads_to(n)))
				if n.mean == None:
					n.mean = 0

			if ancestor_meaning == None:
					self.events['unclear'].append(self.string_branch(self.make_branch_that_leads_to(n)))

			for child in n._children:
				self.events_on_branches(child)

	def event_nodes_maker(self, n):
		if n.is_root() == True:
			self.event_nodes = {'deletion': list(), 'insertion': list(), 'unclear': list()}
			for child in self._children:
				self.event_nodes_maker(child)
		else:
			ancestor_meaning = n._up.mean
			if ancestor_meaning == 1:
				if n.mean == 0:
					self.event_nodes['deletion'].append(n)
				if n.mean == None:
					n.mean = 1

			if ancestor_meaning == 0:
				if n.mean == 1:
					self.event_nodes['insertion'].append(n)
				if n.mean == None:
					n.mean = 0

			if ancestor_meaning == None:
					self.event_nodes['unclear'].append(n)

			for child in n._children:
				self.event_nodes_maker(child)
				
	def events_on_branches_contexts(self, n):
		if n.is_root() == True:
			self.events_contexts = {'deletion': list(), 'insertion': list(), 'unclear': list()}
			for child in self._children:
				self.events_on_branches_contexts(child)
		else:
			ancestor_meaning = n._up.mean
			if ancestor_meaning == 1:
				if n.mean == 0:
					self.events_contexts['deletion'].append(n._up.context)
				if n.mean == None:
					n.mean = 1

			if ancestor_meaning == 0:
				if n.mean == 1:
					self.events_contexts['insertion'].append(n._up.context)
				if n.mean == None:
					n.mean = 0

			if ancestor_meaning == None:
					self.events_contexts['unclear'].append(n._up.context)

			for child in n._children:
				self.events_on_branches_contexts(child)

	def events_on_branches_distances(self, n):
		if n.is_root() == True:
			self.event_distances = {'deletion': list(), 'insertion': list(), 'unclear': list()}
			for child in self._children:
				self.events_on_branches_distances(child)
		else:
			ancestor_meaning = n._up.mean
			if ancestor_meaning == 1:
				if n.mean == 0:
					r = n 
					while r.is_root() != True:
						r = r._up
					d = r.get_distance(n) - (r.get_distance(n, target2=n._up))/2
					self.event_distances['deletion'].append(d)
				if n.mean == None:
					n.mean = 1

			if ancestor_meaning == 0:
				if n.mean == 1:
					r = n 
					while r.is_root() != True:
						r = r._up
					d = r.get_distance(n) - (r.get_distance(n, target2=n._up))/2
					self.event_distances['insertion'].append(d)
				if n.mean == None:
					n.mean = 0

			# if ancestor_meaning == None:
					# self.event_distances['unclear'].append(self.string_branch(self.make_branch_that_leads_to(n)))

			for child in n._children:
				self.events_on_branches_distances(child)


def turn_alignment_to_matrix(in_bp_dict):
	species_list = list(in_bp_dict.keys())

	bp_dict = dict()

	#delete paralouges - we just don't put them in the final martix
	for elem in species_list:
		sp = elem[:4]
		if sp in bp_dict:
			del bp_dict[sp]
		else:
			bp_dict[sp] = in_bp_dict[elem]

	species_list = list(bp_dict.keys())
	species_list.sort()
	matrix = list()

	for species in species_list:
		matrix.append(list(bp_dict[species]))
	matrix = list(map(list, zip(*matrix))) #transposition - now column of the alignment corresponds to a line in the matrix
	return bp_dict, matrix, species_list


def cut_context(pos, species, two_side_length, bp_dict):
	return bp_dict[species][pos - two_side_length:pos + two_side_length+1]


# local_tree = Ancestor_Tree('(petz:0.14339273,(raik:0.26772623,((octo:0.34300359,harp:0.21398805)12:0.08092917,(eury:0.21540617,(rari:0.27565106,(foca:0.26750949,(minu:0.12310897,cras:0.16946154)30:0.12093409)24:0.05101473)9:0.04135668)13:0.03839098)11:0.04916400):0.14339273);')

# print(local_tree.get_tree_root())
# a = 'AAAAA'
# print(local_tree.get_leaves_by_name('foca'))
# print(local_tree.get_leaves_by_name('foca')[0])
# local_tree.get_leaves_by_name('foca')[0].leaf_context('AAAAA')
# local_tree.get_leaves_by_name('minu')[0].leaf_context('AATAA')
# local_tree.get_leaves_by_name('cras')[0].leaf_context('AAAAA')
# local_tree.get_leaves_by_name('harp')[0].leaf_context('AA-AA')
                                                    
# local_tree.node_context_upgoing(5)
# print(local_tree.get_leaves_by_name('cras')[0]._up.context)
# local_tree.node_context_downgoing(5)
# print(local_tree.context)
# print(local_tree.get_leaves_by_name('foca')[0]._up.context)
# print(local_tree.get_leaves_by_name('cras')[0]._up.context)
# print(local_tree.get_leaves_by_name('raik')[0].context)
# print(local_tree.get_leaves_by_name('raik')[0]._up.context)

# print(local_tree.get_distance(local_tree.get_leaves_by_name('cras')[0], 
# 	target2=None, topology_only=False))
