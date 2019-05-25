import numpy as np
from scipy.linalg import fractional_matrix_power as frac_power
eps = pow(10.0, -7)

# Add new point to the end, see how it effects the Circumradius_2^p-Circumradius_1^p quantity.
# I want to test whether [-1, 0, 100] has a larger p_circum_diff or [0, 100]
def p_circumradius_change(alpha,l):
	l2 = list(np.sort(l))
	l1 = l2[:-1]
	circum_2_p = rad_alpha(l2, alpha)
	circum_1_p = rad_alpha(l1, alpha)
	print(circum_2_p)
	print(circum_1_p)
	print(circum_2_p-circum_1_p)
	print("")
	return circum_2_p - circum_1_p
	
def square_circumradius_change(alpha, l):
	l2 = list(np.sort(l))
	l1 = l2[:-1]
	circum_2_p = rad_squared(l2, alpha)
	circum_1_p = rad_squared(l1, alpha)
	print(circum_2_p)
	print(circum_1_p)
	print(circum_2_p-circum_1_p)
	print("")
	return circum_2_p - circum_1_p

def hypothesis_test(l1, new_point, alpha=1.0/3, p = 2):
	old_length = l1[-1] - l1[0]
	new_length = new_point - l1[0]
	l2 = l1 + [new_point]
	circum_2_p = 2*rad_p(l2, alpha, p)
	circum_1_p = 2*rad_p(l1, alpha, p)
	circum_val = circum_2_p-circum_1_p
	length_val = pow(pow(new_length, alpha), p) - pow(pow(old_length, alpha), p)
	print(circum_val)
	print(length_val)
	print (circum_val - length_val)/(length_val)
	return circum_val > length_val - eps

def rad_alpha(ls, alpha=1/3):
	return rad_p(ls, alpha, 1/alpha)

def rad_p(ls, alpha=1/3, p = 3):
	return pow(rad_squared(ls, alpha), p/2.0)

def rad_squared(ls, alpha=1.0/3):
	k = np.linalg.inv(cm(ls, 2*alpha))
	return -k[0][0]/2

def inv_cm(ls, power=2.0/3):
	return np.linalg.inv(cm(ls, power))
#Cayley Menger matrix
def cm(ls, power=2.0/3):
	n = len(ls)+1
	M = np.zeros((n,n))
	D = screw(ls, power)
	for i in range(1, n):
		for j in range(1, n):
			M[i][j] = D[i-1][j-1]
	for i in range(1,n):
		M[0][i]=1
		M[i][0]=1
	return M

def mean_center(n):
	M = np.zeros((n,n))
	for i in range(n):
		for j in range(n):
			M[i][j] = -1.0/n
	for i in range(n):
		M[i][i] += 1
	return M

def gram(ls, power=2.0/3):
	D = screw(ls, power)
	n = len(D)
	M = mean_center(n)
	gram = -np.dot(np.dot(M, D), M)
	return gram 

def schur(M, r):
	B = Bmatrix(M, r)
	Q = res(M, r)
	return Q - np.dot(B, B.T)/M[r][r]
def Bmatrix(M, r):
	n = len(M)
	B = np.zeros((n-1, 1))
	for i in range(n-1):
		if i < r:
			B[i][0] = M[i][r]
		elif i >= r:
			B[i][0] = M[i+1][r]
	return B
def res(M, r):
	n = len(M)
	Mhat = np.zeros((n-1, n-1))
	s = 0; t = 0
	for i in range(n-1):
		for j in range(n-1):
			if i < r:
				s = i
			elif i >= r:
				s = i+1
			if j < r:
				t = j
			elif j >= r:
				t = j+1
			Mhat[i][j] = M[s][t]
	return Mhat

# Distance matrix squared, of screw with screw func |i-j|^power.
def screw(ls, power=2.0/3):
	n = len(ls)
	M = np.zeros((n, n))
	for i in range(n):
		for j in range(n):
			D = abs(ls[i]-ls[j])
			M[i][j] = pow(D, power)
	return M
