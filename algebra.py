import copy


class Element:
	def __init__(self, parent, data):
		self.parent = parent
		self.data = data

	def __str__(self):
		return self.parent.eltstr(self)

	def __eq__(self, other):
		return type(self.parent) == type(other.parent) and self.data == other.data
		# May not be accurate if the same element can be represented different ways


class Module_element(Element):
	def __add__(self, other):
		return self.parent.add(self, other)

	def __rmul__(self, other):
		if type(other) == Algebra_element:
			return self.parent.leftAction(other, self)


class Vector_space_element(Module_element):
	pass


class Algebra_element(Vector_space_element):
	def __mul__(self, other):
		if isinstance(other, Algebra_element) and type(self.parent) == type(other.parent):
			return self.parent.mul(self, other)
		elif isinstance(other, Module_element):
			return other.__rmul__(self)



class GradedF2Module:

	def __init__(self, *args, **kwargs):
		self.basis_elements = []


	# Return a list of basis elements in degree m
	def constructBasis(self, m):
		raise NotImplementedError("Method 'constructBasis' must be implemented in subclasses of GradedF2Module")


	def basisElementPrintableName(self, b):
		raise NotImplementedError("Method 'basisElementPrintableName' must be implemented in subclasses of GradedF2Module")


	def basis(self, m):
		if len(self.basis_elements) <= m:
			self.basis_elements += [None for i in range(m + 1 - len(self.basis_elements))]
		if self.basis_elements[m] == None:
			self.basis_elements[m] = self.constructBasis(m)
		return self.basis_elements[m]


	def element(self, data):
		data = copy.deepcopy(data)

		# Combine terms of the same degree
		data.sort()
		i = 0
		while i < len(data) - 1:
			if data[i][0] == data[i+1][0]:
				for j in range(len(data[i][1])):
					data[i][1][j] = (data[i][1][j] + data[i+1][1][j]) % 2
				data.pop(i+1)
			else:
				i += 1

		# Remove zero terms
		for i in range(len(data) - 1, -1, -1):
			if sum(data[i][1]) == 0:
				data.pop(i)

		return Vector_space_element(self, data)


	def add(self, a, b):
		return self.element(a.data + b.data)


	def eltstr(self, elt):
		s = " + ".join([" + ".join([
			self.basisElementPrintableName(self.basis(term[0])[i]) for i in range(len(self.basis(term[0]))) if term[1][i] == 1
		]) for term in elt.data])
		return s if len(s) > 0 else "0"
