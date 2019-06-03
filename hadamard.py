import numpy as np
def structured_hadamard(d):
	if (d == 0):
		return np.array([[1]])
	A = structured_hadamard(d-1)
	B1 = np.zeros((len(A), 2*len(A)))
	B2 = np.zeros((len(A), 2*len(A)))
	for i in range(len(A)):
		B1[i] = np.concatenate((A[i], A[i]))
		B2[i] = np.concatenate((A[i], -A[i]))
	return np.concatenate((B1, B2))
