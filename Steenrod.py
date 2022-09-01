from algebra import *
import math
import copy

class Steenrod_algebra: # at the prime 2

	def __init__(self, p):
		self.p = 2

	def adem(self, I):
		return Algebra_element(self, self.simplify(copy.deepcopy(I)))

	def simplify(self, I):
		# Remove zeros
		for i in range(len(I)):
			while len(I[i]) > 1 and 0 in I[i]:
				I[i].remove(0)

		# Apply Adem relations
		i = 0
		while i < len(I):
			usedAdem = False
			for j in range(len(I[i]) - 1):
				if I[i][j] < 2 * I[i][j+1]:
					term = I.pop(i)
					for squares in self.adem_relation(term[j], term[j+1]):
						I.append(term[:j] + squares + term[j+2:])
					usedAdem = True
					break
			if not usedAdem: i += 1

		# Remove duplicate terms
		I.sort()
		i = 0
		while i < len(I) - 1:
			if I[i] == I[i+1]:
				I.pop(i)
				I.pop(i)
			else: i += 1

		return I

	def adem_relation(self, i, j):
		removeZeros = lambda X : [x for x in X if x != 0]
		return [removeZeros([i + j - k, k]) for k in range(i // 2 + 1) if combMod2(j - k - 1, i - 2*k) == 1]

	def degree(self, a):
		return [sum(term) for term in a.data]

	def add(self, a, b):
		return self.adem(a.data + b.data)

	def mul(self, a, b):
		return self.adem([c + d for c in a.data for d in b.data])

	def eltstr(self, elt):
		if len(elt.data) == 0: return "0"
		return " + ".join([ " ".join(["Sq^{}".format(i) for i in mon]) for mon in elt.data ])


# efficiently computes nCk for non-negative integers n, k
def combMod2(n, k):
	if k > n: return 0
	n, k = "{0:b}".format(n), "{0:b}".format(k)
	k = "0" * (len(n) - len(k)) + k
	return 1 if all([n[i] == "1" or k[i] == "0" for i in range(len(n))]) else 0



A = Steenrod_algebra(2)
