import copy

class Matrix:
	def __init__(self, ent, mod=False):
		self.m = len(ent)
		self.n = len(ent[0])
		self.ent = copy.deepcopy(ent) # values
		self.mod = mod


	def inv_mod2(self):
		if self.m != self.n:
			raise ValueError("Cannot invert a non-square matrix.")
		n = self.n
		A = [self.ent[r] + [1 if i == r else 0 for i in range(n)] for r in range(n)]
		A = [[A[r][c] % 2 for c in range(2*n)] for r in range(n)]

		for i in range(n):
			if A[i][i] == 0:
				r1 = min([r for r in range(i, n) if A[r][i] == 1])
				for c in range(2*n): A[i][c] = (A[i][c] + A[r1][c]) % 2
			for r in range(i+1, n):
				if A[r][i] == 1:
					for c in range(2*n): A[r][c] = (A[r][c] + A[i][c]) % 2

		for i in range(n-1, 0, -1):
			for r in range(i):
				if A[r][i] == 1:
					for c in range(2*n): A[r][c] = (A[r][c] + A[i][c]) % 2

		return Matrix([[A[r][n + c] for c in range(n)] for r in range(n)], mod=2)


	def transpose(self):
		return Matrix([[self.ent[c][r] for c in range(self.m)] for r in range(self.n)], self.mod)


	def det(self):
		if self.m != self.n:
			raise ValueError("Cannot take the determinant of a non-square matrix.")

		if self.m == 1: return self.ent[0][0]
		ans = sum([self.ent[0][i] * ((-1)**i) * Matrix([[self.ent[j][k] for k in list(range(0, i)) + list(range(i+1, self.m))] for j in range(1, self.m)]).det() for i in range(self.m)])
		return ans if self.mod == False else ans % self.mod


	def __add__(self, other):
		if self.m != other.m or self.n != other.n:
			raise ValueError("Cannot add a {}x{} matrix and a {}x{} matrix.".format(self.m, self.n, other.m, other.n))
		if self.mod != other.mod:
			raise ValueError("Cannot add matrices with entries in different rings.")

		ans = Matrix([[self.ent[r][c] + other.ent[r][c] for c in range(self.n)] for r in range(self.m)], mod=self.mod)
		if self.mod != False:
			for i in range(ans.m):
				for j in range(ans.n):
					ans.ent[i][j] %= self.mod
		return ans


	def __mul__(self, other):
		if self.n != other.m:
			raise ValueError("Cannot multiply a {}x{} matrix and a {}x{} matrix.".format(self.m, self.n, other.m, other.n))
		if self.mod != other.mod:
			raise ValueError("Cannot add matrices with entries in different rings.")

		ans = Matrix([[sum([self.ent[r][i] * other.ent[i][c] for i in range(self.n)]) for c in range(other.n)] for r in range(self.m)], mod=self.mod)
		if self.mod != False:
			for i in range(ans.m):
				for j in range(ans.n):
					ans.ent[i][j] %= self.mod
		return ans


	def __str__(self):
		l = max([len(str(self.ent[r][c])) for c in range(self.n) for r in range(self.m)])
		return "\n".join(["".join([str(self.ent[r][c]) + " "*(l + 1 - len(str(self.ent[r][c]))) for c in range(self.n)]) for r in range(self.m)])
