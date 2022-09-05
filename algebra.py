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
