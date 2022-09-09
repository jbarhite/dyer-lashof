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
