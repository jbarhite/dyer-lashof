from algebra import *
from matrix import *
import Steenrod
import partitions
import copy

class Brown_Gitler_polynomial_algebra:

	def __init__(self):
		self.A = Steenrod.Steenrod_algebra(2)

	# The first entry represents the coefficient, and subsequent entries represent powers of the variables.
	# Example: 1 + x0 x1 ^ 2 + x5  would be input as [[1], [1, 1, 2]], [1, 0, 0, 0, 0, 0, 5]]
	def element(self, data):
		data = [term[0:max([i for i in range(len(term)) if term[i] != 0])+1] for term in copy.deepcopy(data) if term[0] != 0]
		data.sort()
		i = 0
		while i < len(data) - 1:
			if data[i] == data[i+1]:
				data.pop(i)
				data.pop(i)
			else: i += 1
		return Algebra_element(self, data)


	def add(self, a, b):
		return self.element(a.data + b.data)


	def mul(self, a, b):
		product = []
		for t1 in a.data:
			for t2 in b.data:
				product.append([t1[0] * t2[0]] + [(t1[i] if i < len(t1) else 0) + (t2[i] if i < len(t2) else 0) for i in range(1, max(len(t1), len(t2)))])
		return self.element(product)


	def eltstr(self, elt):
		terms = []
		for term in elt.data:
			if term == [1]: terms.append("1")
			else:
				terms.append(" ".join(["x{}^{}".format(i - 1, term[i]) if term[i] > 1 else "x{}".format(i - 1) for i in range(1, len(term)) if term[i] != 0]))
		return "0" if len(terms) == 0 else " + ".join(terms)


	basis_elements = []
	def basis(self, n, m=-1):
		B = Brown_Gitler_polynomial_algebra.basis_elements
		if len(B) < n + 1: B += [False for i in range(n + 1 - len(B))]
		if B[n] == False:
			B[n] = [self.element([[1] + v]) for v in partitions.partitionsByPowersOf2(n, True)]
		return B[n] if m == -1 else [b for b in B[n] if sum(b.data[0][1:]) == m]


	# Sq^a x_i^b = (b choose a) x_i^(b-a) x_(i-1)^(2a), where (b choose a) = 0 if a > b
	def leftAction(self, a, x):
		ans = self.element([])
		for a1 in a.data:
			if a1 == [0]:
				ans += x
				continue
			for x1 in x.data:
				ans1 = [x1]
				for i in reversed(a1):
					ans2 = []
					for term in ans1:
						# Apply Sq^i to term
						if i > sum(term[1:]): continue # instability condition

						xs = [j for j in range(1, len(term)) if term[j] != 0]
						for part in self.parts(i, len(xs)):
							if all([Steenrod.combMod2(term[xs[j]], part[j]) == 1 for j in range(len(xs))]) and not (xs[0] == 1 and part[0] > 0):
								newterm = [1] + [0 for j in range(len(term))]
								for j in range(len(xs)):
									newterm[xs[j]] += term[xs[j]] - part[j]
									newterm[xs[j]-1] += 2 * part[j]
								ans2.append(newterm)
					ans1 = ans2
				ans += self.element(ans1)

		return ans


	conversion_matrices = []

	def elementFromJ(self, x):
		J = x.parent
		d = sum([len(b) for b in J.basis])
		v = Matrix([[0 for i in range(d)]], mod=2)
		for term in x.data:
			w = Matrix([[0 for i in range(sum([len(J.basis[m]) for m in range(term[0])]))] + term[1] + [0 for i in range(sum([len(J.basis[m]) for m in range(term[0] + 1, d + 1)]))]], mod=2)
			v += w

		cm = type(self).conversion_matrices
		if len(cm) < J.n + 1:
			cm += [None for i in range(J.n + 1 - len(cm))]
		if cm[J.n] == None:
			# Produce conversion matrix
			cm[J.n] = Matrix([[0 for c in range(d)] for r in range(d)], mod=2)

		w = cm[J.n] * v.transpose()
		return sum([self.basis(J.n)[i] for i in range(d) if w.v[i] == 1], self.element([]))


	# used in leftAction
	# returns ordered partitions of n into k summands
	def parts(self, n, k):
		if k == 1: return [[n]]
		parts = []
		for i in range(n+1):
			parts += [[i] + part for part in self.parts(n - i, k - 1)]
		return parts


	def printActions(self, n):
		printed = False
		for m in range(min([i for i in range(n + 1) if len(self.basis(i)) > 0]), n + 1):
			if printed and m < n - 1: print()
			for b in self.basis(n, m):
				i = 1
				while i <= m and i + m <= n:
					a = self.A.adem([[i]])
					print("{} * ({}) = {}".format(a, b, a * b))
					printed = True
					i *= 2


class Brown_Gitler_module:

	free_modules = [] # Won't need this if we make Free_unstable_module a singleton class
	instances = []

	def __new__(cls, *args):
		if len(cls.instances) < args[0] + 1:
			cls.instances += [None for i in range(args[0] + 1 - len(cls.instances))]
		if cls.instances[args[0]] == None:
			cls.instances[args[0]] = object.__new__(cls)
			cls.instances[args[0]].init(args[0])
		return cls.instances[args[0]]


	# Gets called directly by __new__
	def init(self, n):
		self.n = n
		fm = Brown_Gitler_module.free_modules
		if len(Brown_Gitler_module.free_modules) < n + 1:
			fm += [Free_unstable_module(i) for i in range(len(fm), n+1)]
		self.basis = [fm[n].basis(m) for m in range(0, n + 1)]
		self.A = Steenrod.Steenrod_algebra(2)


	# returns an object representing an element of J(n)
	# each entry of data is [m, c], where c is a coefficient vector with respect to the basis represented by self.basis(m)

	def element(self, data):
		# Combine terms of the same degree
		data = copy.deepcopy(data)
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
			"({})*".format(self.printFreeModuleBasisElement(self.basis[comp[0]][i])) for i in range(len(self.basis[comp[0]])) if comp[1][i] == 1
		]) for comp in elt.data])
		return s if len(s) > 0 else "0"


	def leftAction(self, a, x):
		ans = self.element([])
		for a1 in a.data:
			for x1 in x.data:
				m, m2 = x1[0], x1[0] + sum(a1)
				if m2 > self.n: continue

				c = [0 for i in range(len(self.basis[m2]))]
				for i in range(len(self.basis[m2])):
					for j in range(len(self.basis[m])):
						if x1[1][j] == 0: continue
						for term in (self.A.adem([self.basis[m2][i][1:]]) * self.A.adem([a1])).data:
							if term == self.basis[m][j][1:]:
								c[i] = (c[i] + 1) % 2
				ans += self.element([[m2, c]])

		return ans


	def mu(self, x, y):
		if x.parent.n + y.parent.n != self.n:
			raise ValueError("Improper inputs to the map J(m) ⊗ J(n) ---> J(m + n).")

		ans = self.element([])
		for t1 in x.data:
			for t2 in y.data:
				c = []
				for b in self.basis[t1[0] + t2[0]]:
					z = [[x, y]]
					for k in b[-1:0:-1]:
						z = [[self.A.adem([[i]]) * term[0], self.A.adem([[k-i]]) * term[1]] for term in z for i in range(k+1)]
						for i in range(len(z) - 1, -1, -1):
							if str(z[i][0]) == "0" or str(z[i][1]) == "0":
								z.pop(i)
					c.append(len(z) % 2)
				ans += self.element([[t1[0] + t2[0], c]])

		return ans


	def elementFromT(self, x):
		if not all([self.n == sum([2**(i-1)*term[i] for i in range(1, len(term))]) for term in x.data]):
			raise ValueError("{} is not an element of J({})".format(x, self.n))

		ans = self.element([])
		for term in x.data:
			ans1 = Brown_Gitler_module(0).element([[0, [1]]])
			m = 0
			for i in range(1, len(term)):
				for j in range(term[i]):
					ans1 = Brown_Gitler_module(m + 2**(i-1)).mu(ans1, Brown_Gitler_module(2**(i-1)).element([[1, [1]]]))
					m += 2**(i-1)
			ans = ans + ans1
		return ans


	def printFreeModuleBasisElement(self, s):
		return "Σ^{} ".format(s[0]) + " ".join(["Sq^{}".format(i) for i in s[1:]])

	
	def printBases(self):
		for m in range(0, self.n + 1):
			if len(self.basis[m]) == 0:
				print("J({})^{} = 0".format(self.n, m))
			else:
				print("J({})^{} = F_2[ {} ]*".format(self.n, m, " , ".join([self.printFreeModuleBasisElement(elt) for elt in self.basis[m]])))


	def printActions(self):
		n = self.n
		for m in range(min([i for i in range(n + 1) if len(self.basis[i]) > 0]), n + 1):
			for b in range(len(self.basis[m])):
				i = 1
				while i <= m and i + m <= n:
					a = self.A.adem([[i]])
					x = self.element([[m, [1 if j == b else 0 for j in range(len(self.basis[m]))]]])
					print("{} * ({}) = {}".format(a, x, a * x))
					i *= 2
			if m < n - 1: print()


class Free_unstable_module:

	instances = []

	def __new__(cls, *args):
		if len(cls.instances) < args[0] + 1:
			cls.instances += [None for i in range(args[0] + 1 - len(cls.instances))]
		if cls.instances[args[0]] == None:
			cls.instances[args[0]] = object.__new__(cls)
			cls.instances[args[0]].init(args[0])
		return cls.instances[args[0]]


	# Gets called directly by __new__
	def init(self, n):
		self.n = n


	def basis(self, m):
		return sorted([[m] + partition for partition in self.parts(self.n - m, self.n // 2)])


	def parts(self, total, maxFirst):
		if 2 * maxFirst < total: return []
	
		ans = []
		for i in range((total + 1) // 2 - 1, min(maxFirst + 1, total + 1)):
			if i == total: ans.append([i])
			else: ans += [[i] + part for part in self.parts(total - i, i // 2)]
		return ans



########################


J = [Brown_Gitler_module(n) for n in range(9)]
T = Brown_Gitler_polynomial_algebra()

x = J[6].element([[2, [1]], [3, [0, 1]]])
print("{} corresponds to {}".format(x, T.elementFromJ(x)))

# y = T.element([[1, 0, 2, 1], [1, 8]])
# print("{} corresponds to {}".format(y, J[8].elementFromT(y)))
