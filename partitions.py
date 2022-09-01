class partitionsByPowersOf2:

	def __init__(self, n):
		if not isinstance(n, int) or n < 1:
			raise ValueError("Input must be a positive integer.")
		self.n = n
		self.p = []


	def __iter__(self):
		return self


	def __next__(self):
		if len(self.p) == 0:
			s = "{0:b}".format(self.n)
			self.p = [2**(len(s) - i - 1) for i in range(len(s)) if s[i] == "1"]
			return self.p

		i = len(self.p) - 1
		while i >= 0 and self.p[i] == 1: i -= 1
		if i == -1:
			raise StopIteration

		self.p = self.p[:i] + [self.p[i] // 2]
		maxpart = self.p[i]
		remaining = self.n - sum(self.p)
		while remaining > 0:
			while maxpart > remaining:
				maxpart //= 2
			self.p.append(maxpart)
			remaining -= maxpart
		return self.p
