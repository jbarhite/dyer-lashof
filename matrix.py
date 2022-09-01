import copy

class Matrix:
	def __init__(self, v):
		self.m = len(v)
		self.n = len(v[0])
		self.v = copy.deepcopy(v) # values


	def inv_mod2(self):
		if self.m != self.n:
			raise ValueError("Cannot invert a non-square matrix.")
		n = self.n
		A = [self.v[r] + [1 if i == r else 0 for i in range(n)] for r in range(n)]
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

		return Matrix([[A[r][n + c] for c in range(n)] for r in range(n)])


	def det(self):
		if self.m != self.n:
			raise ValueError("Cannot take the determinant of a non-square matrix.")
		if self.m == 1: return self.v[0][0]
		return sum([self.v[0][i] * ((-1)**i) * Matrix([[self.v[j][k] for k in list(range(0, i)) + list(range(i+1, self.m))] for j in range(1, self.m)]).det() for i in range(self.m)])


	def __add__(self, other):
		if self.m != other.m or self.n != other.n:
			raise ValueError("Cannot add a {}x{} matrix and a {}x{} matrix.".format(self.m, self.n, other.m, other.n))
		return Matrix([[self.v[r][c] + other.v[r][c] for c in range(self.n)] for r in range(self.m)])


	def __mul__(self, other):
		if self.n != other.m:
			raise ValueError("Cannot multiply a {}x{} matrix and a {}x{} matrix.".format(self.m, self.n, other.m, other.n))
		return Matrix([[sum([self.v[r][i] * other.v[i][c] for i in range(self.n)]) for c in range(other.n)] for r in range(self.m)])


	def __str__(self):
		l = max([len(str(self.v[r][c])) for c in range(self.n) for r in range(self.m)])
		return "\n".join(["".join([str(self.v[r][c]) + " "*(l + 1 - len(str(self.v[r][c]))) for c in range(self.n)]) for r in range(self.m)])
