import numpy as np
from inequality_test import circumcenter_barycentric_coords
# Test whether screw circumcenters are interior
from screw import screw

# I'm testing this since it seems it might not be true.
# It seems true, but note that I haven't unit tested the sub-methods.
def test_screw_circumcenter_interior():
	n = 5 
	p = 2.001
	for i in range(40000):
		ls = np.sort(10*np.random.rand(n))
		D = screw(ls, 2/p)
		b = circumcenter_barycentric_coords(D)
		if(ls[-2] < 0.8 and ls[-1] > 8):
			print(ls)
		if has_neg_coord_barycentric(b):
			print("Distance matrix:")
			print(np.around(D, 3))
			return False
	return True
# To do test: Run test_circumcenter_interior()
# Note that I haven't checked that the overall setup is right, even though the
# submethods seem to be fine.

#This file: wanted to test whether randomly-chosen/structured columns
# of structured hadamard matrices had interior circumcenters. This was since
# the coordinates of the screw simplices, read as columns,
# are column subsets of select hadamard matrices.

#All of my tests, from picking out random hadamard cols
# to picking out structured hadamard
#cols, turn out to be false.

# Now this casts doubt that the circumradius of the screw simplex is always interior.
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

def randomly_weight_rows(M):
	S = M.copy()
	for i in range(len(M)):
		S[i] = np.random.rand()*M[i]
	return S
# Pick out columns 111111, 1111-1-1-1-1, 11-1-111-1-1, 1-11-1...1-1
# Just hacked together for d =3
def pick_out_structured_columns(M):
	n = len(M)
	c = 4
	selected_cols = np.zeros((n,c))
	col_indices = [0, 1, 2, 4]
	for i in range(c):
		ci = col_indices[i]
		selected_cols[:, i] = M[:, ci]
	return selected_cols

def pick_out_random_columns(M, c):
	n = len(M)
	# Picks out a random subset of  size c from np.arange(n)
	col_indices = np.random.choice(n, c, replace = False)
	selected_cols = np.zeros((n,c))
	for i in range(c):
		ci = col_indices[i]
		selected_cols[:, i] = M[:, ci]
	return selected_cols

# When M is a column of points, and M is not necessarily square.
def dist_mat_squared(M):
	n = len(M[0])
	A = np.zeros((n, n))
	for i in range(n):
		for j in range(n):
			A[i][j] = np.linalg.norm(M[:,i]-M[:,j]) ** 2
	return A

def test_circumcenter_interior():
	n = 3
	c = 4
	for i in range(1000):
		M = structured_hadamard(n)
		structured_hadamard_cols = pick_out_structured_columns(M)
		P = randomly_weight_rows(structured_hadamard_cols)
		#randomly_chosen_hadamard_columns = pick_out_random_columns(structured_hadamard(n), c)
		#P = randomly_weight_rows(randomly_chosen_hadamard_columns)
		D = dist_mat_squared(P)
		b = circumcenter_barycentric_coords(D)
		if has_neg_coord_barycentric(b):
			print("Points are:")
			print(np.around(P, 3))
			print("Distance matrix:")
			print(np.around(D, 3))
			return False
	return True


def has_neg_coord_barycentric(b):
	eps = 0.0001
	if (sum(b) < 1-eps or sum(b) > 1+eps):
		print("Not a barycentric coord. Sum is not right.")
		print(sum(b))
		return True
	for e in b:
		if e < - eps:
			print("Barycentric coord has a negative coeff")
			print(e)
			return True
	return False

