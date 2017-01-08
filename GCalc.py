DEFAULT_IRR = 0x1b
MIN_N = 2
DEFAULT_N = 8
MAX_N = 11

class GCalc:

	def __init__(self, N=8, irr=0x1b, check=True):
		if check and (N < MIN_N or N > MAX_N):
			print("N must be a number in the range " + str(MIN_N) + " - " + str(MAX_N) + " inclusive, setting to default value: " + str(DEFAULT_N))
			N = DEFAULT_N
		self.N = N
		self.pN = 2**N
		self.poly = irr & (self.pN - 1)

	def getN(self):
		return self.N

	def getpN(self):
		return self.pN

	def getIrreducible(self):
		return self.poly

	def add(self, x, y):
		if (x < 0 or x >= self.pN or y < 0 or y >= self.pN):
			print("One or more supplied values was not an integer in allowed range: 0-" + str(self.pN - 1))
			return -1

		return (x ^ y)

	def gmult2(self, x):
		if (x < 0 or x >= self.pN):
			print("Supplied value was not an integer in the allowed range 0-" + str(self.pN - 1))
			return -1

		return ((x << 1) ^ ( ( (x >> (self.N - 1) ) & 0x01) * self.poly)) & (self.pN -1)

	def mult(self, x, y):
		if (x < 0 or x >= self.pN or y < 0 or y >= self.pN):
			print("One or more supplied values was not an integer in allowed range 0-" + str(self.pN - 1))
			return -1

		m = (y & 0x01) * x 
		for bit in range(1, self.N):
			m2 = self.gmult2(x)
			for shift in range(1, bit):
				m2 = self.gmult2(m2)

			m ^= ((y >> bit) & 0x01 ) * m2

		return m & (self.pN - 1)

	def pow(self, x, y):
		if (x < 0 or x >= self.pN):
			print("Supplied value outside allowed range 0-" + str(self.pN - 1))
			return -1

		if y == 0:
			return 1

		power = x
		if y > 0:
			for i in range(1, y):
				power = self.mult(power, x)
			
		if y < 0:
			count = 0 - y
			for i in range(1, count):
				power = self.mult(power, x)
			power = self.pow(power, self.pN - 2)

		return power

	def div(self, x, y):
		if (x < 0 or x >= self.pN or y < 0 or y >= self.pN):
			print("One or more supplied values was not an integer in allowed range 0-" + str(self.pN - 1))
			return -1

		yinverse = self.pow(y, -1)
		if (yinverse < 0):
			return -1
		return self.mult(x, yinverse)
