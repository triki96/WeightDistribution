# Read a matrix from file
def readFromFile(filename):
	with open(filename) as f:
		array = [[int(x) for x in line.split()] for line in f]
	# Read code parameters
	n = array[0][0] # number of cols of the parity check matrix
	m = array[0][1] # number of rows of the parity check matrix
	k = n-m
	# Position 1 of array denotes the max weight of rows and columns
	# Position 2 of array denotes the weight of every rows and columns
	# Positions 4 to n+4 denotes the 1's position in every rows
	H = zero_matrix(GF(2),m, n)
	for j in range(n):
		for i in range(len(array[4+j])):
			H[array[j+4][i]-1,j] = 1
	return n,m,H
