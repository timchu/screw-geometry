import numpy as np

def laplacian(A):
	n = len(A)
	D = np.zeros((n,n))
	for i in range(n):
		rowsum = 0 
		for j in range(n):
			rowsum += A[i][j]
		D[i][i] = rowsum
	return D-A

def weighted_path_adjacency(n):
	A = np.zeros((n,n))
	for i in range(1, n):
		w = 2*np.random.rand()
		A[i][i-1]=w
		A[i-1][i]=w
	return A

def path_adjacency(n):
	A = np.zeros((n,n))
	for i in range(1, n):
		A[i][i-1]=1
		A[i-1][i]=1
	return A

def permuted_cycle_adjacency(n):
	ls = np.random.permutation(n)
	A = np.zeros((n,n))
	for q in range(1, n):
		i = ls[q]
		i0 = ls[q-1]
		A[i][i0]=1
		A[i0][i]=1
	i = ls[n-1]
	i0 = ls[0]
	A[i][i0]=1
	A[i0][i]=1
	return A

def cycle_adjacency(n):
	A = np.zeros((n,n))
	for i in range(1, n):
		A[i][i-1]=1
		A[i-1][i]=1
	A[n-1][0]=1
	A[0][n-1]=1
	return A

def random_graph_adjacency(n, m):
	A = np.zeros((n,n))
	for _ in range(m):
		(a, b) = (np.random.randint(0, n), np.random.randint(0, n))
		while a==b:
			(a, b) = (np.random.randint(0, n), np.random.randint(0, n))
		A[a][b] = 1
		A[b][a] = 1
	# The cycle forces connectivity.
	return A+cycle_adjacency(n)

def effective_resistance(A, e1, e2):
	n = len(A)
	Lpinv = np.linalg.pinv(laplacian(A))
	chi = np.zeros(n)
	chi[e1] = 1
	chi[e2]=-1
	return chi @ Lpinv @ chi
def neg_er(A, e1, e2):
	return -effective_resistance(A, e1, e2)

def effective_conductance(A, e1, e2):
	return 1/effective_resistance(A, e1, e2)

def resistance_distances(A):
	n = len(A)
	lls = [[effective_resistance(A, i, j) for i in range(n)] for j in range(n)]
	return np.array(lls)

def inverse_ones_times_resistance(A):
	n = len(A)
	ones = np.ones(n)
	return nk

# The input Func takes in adjacency matrix and e1, e2
def concavity_test(func):
	for i in range(1000):
		print(i)
		n = 10
		A = 2*path_adjacency(n)
		B = random_graph_adjacency(n, n)
		(A,B)=(permuted_cycle_adjacency(n), permuted_cycle_adjacency(n))
		(e1, e2) = (0, n-1)
		# If: convex anywhere, false.
		convex_val = func(A, e1, e2)+func(B, e1, e2) - 2*func((A+B)/2, e1, e2)
		print(convex_val)
		print (func(A, e1, e2))
		print (func(B, e1, e2))
		print((func(A+B)/2, e1, e2))
		if convex_val > 0.00001:
			return False
	return True
# This should return true
def neg_ER_concavity_test():
	return concavity_test(neg_er)
def ER_concavity_test():
	return concavity_test(effective_resistance)
def EC_concavity_test():
	return concavity_test(effective_conductance)


def adjacency_to_list(A):
	return np.hstack(A)

# Requires the length of ls be a square
def list_to_adjacency(ls):
	n = int(np.sqrt(len(ls)+0.001))
	A = np.zeros((n,n))
	c = 0
	for i in range(n):
		for j in range(n):
			A[i][j] = ls[c]
			c += 1
	return A

def from_list_effective_resistance(ls, e1, e2):
	A = list_to_adjacency(ls)
	return effective_resistance(A, e1, e2)

def from_list_effective_conductance(ls, e1, e2):
	A = list_to_adjacency(ls)
	return effective_conductance(A, e1, e2)