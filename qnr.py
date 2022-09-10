import Steenrod
import Brown_Gitler as BG


class Polynomial_Ring_A_Module(Steenrod.Graded_A_module):

	def constructBasis(self, m):
		return [m] if m >= 0 else []


	def basisElementPrintableName(self, b):
		if b == 0: return "1"
		if b == 1: return "t"
		return "t^{}".format(b)


	def leftActionOnBasis(self, i, m, j):
		return self.element([[m + i, [Steenrod.combMod2(m, i)]]])


class Polynomial_Ring_A_Module_t_Inverse(Steenrod.Graded_A_module):

	def constructBasis(self, m):
		return [m] if m >= -1 else []


	def basisElementPrintableName(self, b):
		if b == -1: return "t^-1"
		if b == 0: return "1"
		if b == 1: return "t"
		return "t^{}".format(b)


	def leftActionOnBasis(self, i, m, j):
		if m == -1: return self.element([[m + i, [1]]])
		return self.element([[m + i, [Steenrod.combMod2(m, i)]]])


n, r = 2, 4
J = [BG.Brown_Gitler_module(i) for i in range(max(n, n + r) + 1)]
F2t = Polynomial_Ring_A_Module()
M = Steenrod.Graded_A_tensor_product(J[n], F2t)
T = BG.Brown_Gitler_polynomial_algebra()

def sumOfBasisElements(J, M, n, r, z):
	if len(z.data) == 0: return 0
	ans = 0
	for i in range(len(z.data[0][1])):
		if z.data[0][1][i] == 1:
			a = M.basis(n + r)[i][0]
			ans += len(T.elementFromJ(J.basisElement(a[0], a[1])).data)
	return ans % 2

f = lambda z : sumOfBasisElements(J[n], M, n, r, z)

for m in range(n + 1):
	for b in range(len(T.basis(n, m))):
		for t in range(n + r - m + 1):
			x = M.tensor(J[n].elementFromT(T.basis(n, m)[b]), F2t.basisElement(t, 0))
			print("{} ⊗ {} maps to {}".format(
				T.basis(n, m)[b],
				F2t.basisElement(t, 0),
				T.elementFromJ(J[n + r].mapIntoJ(f, x, m + t))
			))
