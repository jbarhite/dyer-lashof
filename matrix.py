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


	def row_reduce_mod2(self):
		m, n = self.m, self.n

		# first make sure that every entry is reduced mod 2
		for r in range(m):
			for c in range(n):
				self.ent[r][c] = self.ent[r][c] % 2

		pivotRow = 0 # this is the row in which we are currently trying to produce a pivot
		pivotCols = []
		for c in range(self.n):
			for r in range(pivotRow, self.m):
				# find the first nonzero entry below (or in) the pivot position
				# if none of these entries are nonzero, then nothing needs to be done to this column
				if self.ent[r][c] != 0:
					# swap rows to bring a nonzero entry (i.e. 1) into the pivot position
					self.swapRows(pivotRow, r)

					# clear entries below the pivot
					for r2 in range(pivotRow + 1, self.m):
						if self.ent[r2][c] != 0:
							self.addMultOfRow(r2, pivotRow, 1, mod2=True)
					
					pivotRow += 1
					pivotCols.append(c)
					break

		# the matrix is now in row echelon form

		for r in reversed(range(len(pivotCols))):
			for r2 in range(r):
				if self.ent[r2][pivotCols[r]] != 0:
					self.addMultOfRow(r2, r, 1, mod2=True)

		# the matrix is now in reduced row echelon form


	# WARNING: this method assumes that the matrix is already in row echelon form
	def rank(self):
		return max(r for r in range(self.m) if any(self.ent[r][c] != 0 for c in range(self.n))) + 1


	# replaces row r1 with row r1 + k times row r2
	def addMultOfRow(self, r1, r2, k, mod2=False):
		for c in range(self.n):
			newEntry = self.ent[r1][c] + k * self.ent[r2][c]
			if mod2: newEntry = newEntry % 2
			self.ent[r1][c] = newEntry


	def swapRows(self, r1, r2):
		temp = [self.ent[r1][c] for c in range(self.n)]
		for c in range(self.n):
			self.ent[r1][c] = self.ent[r2][c]
			self.ent[r2][c] = temp[c]


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
