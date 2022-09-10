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



class Graded_F2_module:

	def __init__(self, *args, **kwargs):
		self.basis_elements = []
		self.basis_elements_neg = []


	# Return a list of basis elements in degree m
	def constructBasis(self, m):
		raise NotImplementedError("Method 'constructBasis' must be implemented in subclasses of Graded_F2_module")


	def basisElementPrintableName(self, b):
		raise NotImplementedError("Method 'basisElementPrintableName' must be implemented in subclasses of Graded_F2_module")


	def basis(self, m):
		if m >= 0:
			if len(self.basis_elements) <= m:
				self.basis_elements += [None for i in range(m + 1 - len(self.basis_elements))]
			if self.basis_elements[m] == None:
				self.basis_elements[m] = self.constructBasis(m)
			return self.basis_elements[m]
		else:
			if len(self.basis_elements_neg) < -m:
				self.basis_elements_neg += [None for i in range(-m - len(self.basis_elements_neg))]
			if self.basis_elements_neg[-m - 1] == None:
				self.basis_elements_neg[-m - 1] = self.constructBasis(m)
			return self.basis_elements_neg[-m - 1]


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


	def basisElement(self, m, j):
		return self.element([[m, [1 if i == j else 0 for i in range(len(self.basis(m)))]]])


	def add(self, a, b):
		return self.element(a.data + b.data)


	def eltstr(self, elt):
		s = " + ".join([" + ".join([
			self.basisElementPrintableName(self.basis(term[0])[i]) for i in range(len(self.basis(term[0]))) if term[1][i] == 1
		]) for term in elt.data])
		return s if len(s) > 0 else "0"


	def printBases(self, minIndex, maxIndex):
		for m in range(minIndex, maxIndex + 1):
			if len(self.basis(m)) == 0:
				print("Degree {}: None".format(m))
			else:
				print("Degree {}: {}".format(m, ", ".join([str(self.basisElement(m, i)) for i in range(len(self.basis(m)))])))


class Graded_F2_tensor_product(Graded_F2_module):

	def __init__(self, M, N, minIndex=0):
		self.M = M
		self.N = N
		self.minIndex = minIndex
		super().__init__()


	def constructBasis(self, m):
		return sum([[((k, i), (m - k, j)) for i in range(len(self.M.basis(k))) for j in range(len(self.N.basis(m - k)))] for k in range(self.minIndex, m + 1 - self.minIndex)], [])
		# return sum([[(a, b) for a in self.M.basis(i) for b in self.N.basis(m - i)] for i in range(self.minIndex, m + 1 - self.minIndex)], [])


	def basisElementPrintableName(self, b):
		return "{} ⊗ {}".format(
			self.M.basisElementPrintableName(self.M.basis(b[0][0])[b[0][1]]),
			self.N.basisElementPrintableName(self.N.basis(b[1][0])[b[1][1]])
		)


	def tensor(self, x, y):
		ans = self.element([])
		for t1 in x.data:
			for t2 in y.data:
				for i in range(len(t1[1])):
					for j in range(len(t2[1])):
						if t1[1][i] == 1 and t2[1][j] == 1:
							ans += self.basisElement(t1[0] + t2[0], self.basis(t1[0] + t2[0]).index(((t1[0], i), (t2[0], j))))
		return ans
