import gcalcn
import copy

class GMatrix:

	DIMENSION_MIN_COLUMN=1
	DIMENSION_MIN_ROW=1
	DIMENSION_MAX=64

	def __init__(self, N=8, poly=0x1b):
		self.calc = gcalcn.GCalcN(N, poly)
		self.pN = self.calc.pN
		self.N = self.calc.N

	def verifyGMatrix(self, matrix):
		length = len(matrix[0])

		for row in matrix:
			if len(row) != length:
				print('badly formed matrix')
				return False

		for row in matrix:
			for e in row:
				if e < 0 or e >= self.pN:
					print("Invalid elements! Elements in matrix must be equal to or greater than 0, and less than " + str(self.pN))
					return False
		return True

	def verifyGVector(self, vector):
		for e in vector:
			if e < 0 or e >= self.pN:
					print("Elements in vector must be equal to or greater than 0, and less than " + str(self.pN))
					return False
		return True

	def genCirculant(self, elementVector):
		if not self.verifyGVector(elementVector):
			print("vector contained invalid elements")
			return -1

		if len(elementVector) < self.DIMENSION_MIN_COLUMN or len(elementVector) > self.DIMENSION_MAX:
			print("element vector must have from 1 to 8 elements")
			return -1

		c = []
		for r in range(0, len(elementVector)):
			c.append( elementVector[len(elementVector) - r: len(elementVector)] + elementVector[0: len(elementVector) - r] )
		return c

	def genZeroMatrix(self, rows, columns):
		if rows < self.DIMENSION_MIN_ROW or rows > self.DIMENSION_MAX:
			print("invalid row value, matrix may be 1x1 -> 8x8")
			return -1

		if columns < self.DIMENSION_MIN_COLUMN or columns > self.DIMENSION_MAX:
			print("invalid column value, matrix may be 1x1 -> 8x8")
			return -1
		return [[0 for i in range(0, columns)] for r in range(0, rows)]

	def genRowEchelon(self, rows, columns):
		m = self.genZeroMatrix(rows, columns)
		for i in range(0, rows):
			m[i][i] = 1
		return m

	def genValueMatrix(self, rows, columns, elementlist):
		if not self.verifyGVector(elementlist):
			print("elementlist contained invalid elements")
			return -1

		m = self.genZeroMatrix(rows, columns)
		row = 0
		index = 0
		length = rows * columns
		if len(elementlist) < length:
			length = len(elementlist)
		if m:
			while index < length:
				m[row][index % columns] = elementlist[index]
				index = index + 1
				if index % columns == 0:
					row = row + 1
		return m

	def matrixAdd(self, matrix1, matrix2):
		if not self.verifyGMatrix(matrix1) or not self.verifyGMatrix(matrix2):
			print("one or more matrices were malformed or contained invalid elements")
			return -1

		if len(matrix1) != len(matrix2) or len(matrix1[0]) != len(matrix2[0]):
			print("matrices must have the same dimension")
			return -1
		result = [[0 for c in range(0, len(matrix1[0]))] for r in range(0, len(matrix1)) ]

		for row in range(0, len(result)):
			for column in range(0, len(result[1])):
				result[row][column] = self.calc.add(matrix1[row][column], matrix2[row][column])

		return result

	def matrixProduct(self, matrix1, matrix2):
		if not self.verifyGMatrix(matrix1) or not self.verifyGMatrix(matrix2):
			print("one or more matrices were malformed or contained invalid elements")
			return -1

		if len(matrix1[0]) != len(matrix2):
			print("matrix1 column length must be equal to number of matrix2 rows")
			return -1

		result = [[0 for c in range(0, len(matrix2[0]))] for r in range(0, len(matrix1))]

		for row in range(0, len(result)):
			for column in range(0, len(result[0])):
				result[row][column] = self.vectorProduct(self.getRowVector(matrix1, row), self.getColumnVector(matrix2, column))

		return result

	def vectorProduct(self, v1, v2):

		if not self.verifyGVector(v1) or not self.verifyGVector(v2):
			print("one or more vectors were malformed or contained invalid elements")
			return -1

		if len(v1) != len(v2):
			print("vectors must be the same length")
			return -1

		p = 0
		for e in range(0, len(v1)):
			p = self.calc.add(p, self.calc.mult(v1[e], v2[e]))
		return p

	def vectorAdd(self, v1, v2):
		if not self.verifyGVector(v1) or not self.verifyGVector(v2):
			print("one or more vectors were malformed or contained invalid elements")
			return -1

		if len(v1) != len(v2):
			print("vectors must be the same length")
			return -1

		vResult = []
		for e in range(0, len(v1)):
			vResult.append(self.calc.add(v1[e], v2[e]))

		return vResult

	def getColumnVector(self, matrix, column):

		if not self.verifyGMatrix(matrix):
			print("matrix was malformed or contained invalid elements")
			return -1

		vector = []
		for row in range(0, len(matrix)):
			vector.append(matrix[row][column])
		return vector

	def getRowVector(self, matrix, row):
		return matrix[row]

	def scalarMult(self, m, matrix):
		mult = [[0 for c in range(0, len(matrix[0]))] for r in range(0, len(matrix))]
		for row in range(0, len(matrix)):
			for column in range(0, len(matrix)):
				mult[row][column] = self.calc.mult(m, matrix[row][column])
		return mult


	def getDetR(self, matrix):  #simple recursion, very slow, don't use on matrices larger than 8x8

		if not self.verifyGMatrix(matrix):
			print("matrices was malformed or contained invalid elements")
			return -1

		if len(matrix) != len(matrix[0]):
			print("Matrix must be square")
			return -1

		if len(matrix) == 1:
			return matrix[0][0]

		if len(matrix) == 2:
			return self.calc.add(self.calc.mult(matrix[0][0], matrix[1][1]), self.calc.mult( matrix[0][1], matrix[1][0]))

		det = 0
		for column in range(0, len(matrix)):
			det = self.calc.add( det, self.calc.mult(matrix[0][column], self.getDet(self.getCofactor(column, 0, matrix))) )

		return det

	def getAdj(self, matrix): 
		if not self.verifyGMatrix(matrix):
			print("matrix was malformed or contained invalid elements")
			return -1

		#pretty sure matrix must be square
		if len(matrix) != len(matrix[0]):
			print("Matrix must be square")
			return -1
		c = self.getCofactorMatrix(matrix)
		adj = self.getTranspose(c)
		return adj    

	def getInverseD(self, matrix):  #I know there are much faster ways, don't use this one on matrices larger than 8x8, still working on row reduction determinant

		if not self.verifyGMatrix(matrix):
			print("matrix was malformed or contained invalid elements")
			return -1

		d = self.getDet(matrix)
		if d != 0 and d != -1:
			c = self.getCofactorMatrix(matrix) #according to my readings, thought it was supposed to be the adjugate matrix, but, in fact, it appears to be the cofactor matrix
			return self.scalarMult(self.calc.pow(d, -1), c)
		else:
			print("Determinant is 0 or not defined, inverse does not exist")
			return -1

	def getTranspose(self, matrix):
		if not self.verifyGMatrix(matrix):
			print("matrix was malformed or contained invalid elements")
			return -1

		t = [[0 for c in range(0, len(matrix))] for r in range(0, len(matrix))]
		for row in range(0, len(matrix)):
			for column in range(0, len(matrix)):
				t[column][row] = matrix[row][column]
		return t

	def getCofactorMatrix(self, matrix):
		if not self.verifyGMatrix(matrix):
			print("matrix was malformed or contained invalid elements")
			return -1

		c = [[0 for c in range(0, len(matrix[0]))] for r in range(0, len(matrix))]
		for row in range(0, len(c)):
			for column in range(0, len(c)):
				c[row][column] = self.getDet(self.getCofactor(row, column, matrix))
		return c

	def getCofactor(self, column, row, matrix):
		if not self.verifyGMatrix(matrix):
			print("matrix was malformed or contained invalid elements")
			return -1
		if len(matrix) != len(matrix[0]):
			print("matrix must be square")
			return -1

		if row >= len(matrix) or row < 0:
			print("row out of bounds")
			return -1
		if column >= len(matrix[0]) or column < 0:
			print("column out of bounds")
			return -1 

		c = []
		for r in range(0, len(matrix)):
			if r != row:
				c.append(matrix[r][0:column] + matrix[r][column+1:len(matrix)])

		return c

	#matrix operations
	def rowValueMult(self, m, row, value):
		matrix = copy.deepcopy(m)
		if not self.verifyGMatrix(matrix):
			print("matrix was malformed or contained invalid elements")
			return -1

		if value < 0 or value >= self.pN:
			print("value out of allowed range")
			return -1

		for e in range(0, len(matrix[row])):
			matrix[row][e] = self.calc.mult(matrix[row][e], value)
		return matrix

	def rowSwap(self, m, row1, row2):
		matrix = copy.deepcopy(m)
		if row1 == row2:
			return matrix
		if not self.verifyGMatrix(matrix):
			print("matrix was malformed or contained invalid elements")
			return -1

		if row1 >= len(matrix) or row1 < 0 or row2 >= len(matrix) or row2 < 0:
			print("one of the selected rows was out of range")
			return -1


		x = matrix[row1]
		matrix[row1] = matrix[row2]
		matrix[row2] = x

		return matrix

	def rowAdd(self, m, targetrow, addrow):
		matrix = copy.deepcopy(m)
		if not self.verifyGMatrix(matrix):
			print("matrix was malformed or contained invalid elements")
			return -1

		if targetrow >= len(matrix) or targetrow < 0 or addrow >= len(matrix) or addrow < 0:
			print("one of the selected rows was out of range")
			return -1

		newRow = self.vectorAdd(matrix[targetrow], matrix[addrow])
		if newRow == -1:
			return -1

		matrix[targetrow] = newRow
		return matrix

	def getRowReduced(self, m):
		matrix = copy.deepcopy(m)
		if not self.verifyGMatrix(matrix):
			print("matrix was malformed or contained invalid elements")
			return -1

		pivotC = 0
		rowT = 0

		while pivotC < len(matrix) and pivotC < len(matrix[0]) and rowT < len(matrix):
			#find next pivot
			pivotRow = rowT
			while matrix[pivotRow][pivotC] == 0:
				pivotRow += 1
				if pivotRow == len(matrix):
					pivotC += 1
					if pivotC == len(matrix[0]): #if you are out of columns, return
						return matrix
					pivotRow = rowT

			matrix = self.rowSwap(matrix, rowT, pivotRow)

			#reduce other values in pivot to zero
			for row in range(0, len(matrix)):
				if row != rowT:
					if matrix[row][pivotC] != 0:
						value = self.calc.div(matrix[row][pivotC], matrix[rowT][pivotC])
						matrix = self.rowValueMult(matrix, rowT, value)
						matrix = self.rowAdd(matrix, row, rowT)

			#set pivot value to 1
			if matrix[rowT][pivotC] != 1:
				matrix = self.rowValueMult(matrix, rowT, self.calc.pow(matrix[rowT][pivotC], -1))

			pivotC += 1
			rowT += 1

		#see if it is in row echelon
		for check in range(0, len(matrix)):
			if matrix[check][check] != 1 or sum(matrix[check][0:len(matrix)]) != 1:
				return matrix, False

		return matrix, True

	def getInverse(self, m):
		matrix = copy.deepcopy(m)
		if not self.verifyGMatrix(matrix):
			print("matrix was malformed or contained invalid elements")
			return -1

		if len(matrix) != len(matrix[0]):
			print("Matrix must be square")
			return -1
		re = self.genRowEchelon(len(matrix), len(matrix))
		inv = []
		for row in range(0, len(matrix)):
			inv.append(matrix[row] + re[row])

		inv, isRowForm = self.getRowReduced(inv)
		if not isRowForm:
			print("matrix has no inverse")
			return -1

		for row in range(0, len(matrix)):
			inv[row] = inv[row][len(matrix):]
		return inv

	def getDet(self, m):
		matrix = copy.deepcopy(m)
		if not self.verifyGMatrix(matrix):
			print("matrix was malformed or contained invalid elements")
			return -1

		if len(matrix) != len(matrix[0]):
			print("Matrix must be square")
			return -1

		if len(matrix) == 1:
			return matrix[0][0]

		if len(matrix) == 2:
			return self.calc.add(self.calc.mult(matrix[0][0], matrix[1][1]), self.calc.mult( matrix[0][1], matrix[1][0]))

		mult = []

		pivotC = 0
		rowT = 0

		while pivotC < len(matrix[0]) and rowT < len(matrix):
			for row in matrix:
				if sum(row) == 0:
					return 0
			#find next pivot
			pivotRow = rowT
			while matrix[pivotRow][pivotC] == 0:
				pivotRow += 1
				if pivotRow == len(matrix):
					return 0

			if (rowT != pivotRow):
				matrix = self.rowSwap(matrix, rowT, pivotRow)
				mult.append(1)

			#reduce other values in pivot to zero
			for row in range(rowT + 1, len(matrix)):
				if row != rowT:
					if matrix[row][pivotC] != 0:
						value = self.calc.div(matrix[row][pivotC], matrix[rowT][pivotC])
						mult.append(value)
						matrix = self.rowValueMult(matrix, rowT, value)
						matrix = self.rowAdd(matrix, row, rowT)

			pivotC += 1
			rowT += 1

		det = 1

		for index in range(0, len(matrix)):
			det = self.calc.mult(det, matrix[index][index])

		b = 1
		for value in mult:
			b = self.calc.mult(b, value)

		det = self.calc.div(det, b)

		return det






