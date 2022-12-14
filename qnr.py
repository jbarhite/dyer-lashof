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


def sumOfBasisElements(J, M, n, r, z):
	if len(z.data) == 0: return 0
	ans = 0
	for i in range(len(z.data[0][1])):
		if z.data[0][1][i] == 1:
			a = M.basis(n + r)[i][0]
			ans += len(T.elementFromJ(J.basisElement(a[0], a[1])).data)
	return ans % 2


def alpha(n): return sum([int(i) for i in str(bin(n))[2:]])


def Q1homologyBasis(n):
	if n < 2: raise ValueError("Q1homologyBasis currently only accepts n >= 2")
	if n % 2 == 0:
		return [
			T.element([[1] + [2 * int(i) for i in str(bin(n // 2))[2:][::-1]]]),
			T.element([[1] + [2 * int(i) for i in str(bin(n // 2 - 1))[2:][::-1]]]) * T.element([[1, 0, 1]])
		], [2 * alpha(n // 2), 2 * alpha(n // 2 - 1) + 1]
	if n % 2 == 1:
		return [
			T.element([[1] + [2 * int(i) for i in str(bin(n // 2))[2:][::-1]]]) * T.element([[1, 1]]),
			T.element([[1] + [2 * int(i) for i in str(bin(n // 2 - 1))[2:][::-1]]]) * T.element([[1, 1, 1]])
		], [2 * alpha(n // 2) + 1, 2 * alpha(n // 2 - 1) + 2]


def qnrIsQ1Iso(n, r):
	f = lambda z : sumOfBasisElements(J[n], M[n], n, r, z)
	Q1hb, deg = Q1homologyBasis(n)
	deg = [d + 2 for d in deg]
	Q1hb2, deg2 = Q1homologyBasis(n + r)
	if min(deg) != min(deg2) or max(deg) != max(deg2): return False

	basis = [M[n].tensor(J[n].elementFromT(b), F2t.basisElement(2, 0)) for b in Q1hb] # basis of H(J(n) ??? (t); Q1)
	imageOfBasis = [J[n+r].mapIntoJ(f, basis[i], deg[i]) for i in range(2)]
	zero = J[n+r].element([])
	return imageOfBasis[0] != zero and imageOfBasis[1] != zero and imageOfBasis[0] != imageOfBasis[1]


def degOfJnQ1Homology(n):
	if n <= 1: return (1,)
	return (2 * alpha(n // 2) + (n % 2), 2 * alpha(n // 2) + 1 + (n % 2))


N, R = 30, 30
J = [BG.Brown_Gitler_module(i) for i in range(N + R + 1)]
F2t = Polynomial_Ring_A_Module()
M = [Steenrod.Graded_A_tensor_product(J[n], F2t) for n in range(N + 1)]
T = BG.Brown_Gitler_polynomial_algebra()


for n in range(20):
	print(n, degOfJnQ1Homology(n))
exit()

for n in range(2, N + 1):
	for r in range(R + 1):
		if r % 2 == 1: continue
		if qnrIsQ1Iso(n, r): print("q({}, {}))_* is an isomorphism".format(n, r))
