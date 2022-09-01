class partitionsByPowersOf2:

	def __init__(self, n, outputCoefficientVector=False):
		if not isinstance(n, int) or n < 0:
			raise ValueError("Input must be a non-negative integer.")
		self.n = n
		self.p = []
		self.oCV = outputCoefficientVector


	def __iter__(self):
		return self


	def __next__(self):
		if len(self.p) == 0:
			self.p = [int(x) for x in "{0:b}".format(self.n)[::-1]]
			return self.p[:] if self.oCV else [2**i for i in range(len(self.p)) if self.p[i] == 1][::-1]

		if self.p[0] == self.n:
			raise StopIteration

		l = len(self.p)
		for i in range(1, l):
			if self.p[i] != 0:
				self.p = [0 for j in range(i)] + ([self.p[i] - 1] + self.p[i+1:] if i < l - 1 or self.p[i] > 0 else [])
				r = self.n - sum([2**j * self.p[j] for j in range(i, len(self.p))])
				for j in range(i - 1, -1, -1):
					while r >= 2**j:
						self.p[j] += 1
						r -= 2**j

				if self.oCV: return self.p[:]
				return sum([[2**i for j in range(self.p[i])] for i in range(len(self.p))], [])[::-1]
