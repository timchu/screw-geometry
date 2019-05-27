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

# Tests whether the circumradius of the average of lists l1 and l2, is larger than the average of the circumradii.
def general_convex_test(l1, l2, p):
	l3 = [(l1[i]+l2[i])/2 for i in range(len(l1))]
	diam_1 = diam_p(l1, p)
	diam_2 = diam_p(l2, p)
	diam_3 = diam_p(l3, p)
	# I'm testing whether val is negative. If it's negative, my hypothesis is true, and if not, my hypothesis is false.
	val = diam_1 + diam_2 - 2*diam_3
	if np.random.rand() > 0.99:
		print(diam_1)
		print(val)
	return val <= 0.000001

def mass_general_convex_test(p):
	mass = 1000
	n = 10
	for i in range(mass):
		l1 = list(np.sort(np.random.rand(n)))
		l2 = list(np.sort(np.random.rand(n)))
		if not(general_convex_test(l1, l2, p)):
			print((l1, l2))
			return False
	for i in range(mass):
		l1 = list(np.linspace(0, n, n))+ [n+n*4*(i+1)/100]
		l2 = list(np.sort(np.random.rand(n+1)))
		if not(general_convex_test(l1, l2, p)):
			print("Broke at: ")
			print((l1, l2))
			return False
	l1 = [1,4,5,6,6.1,6.2, 9]
	l2 = list(np.sort([-3,4,-10,5,6,1,1.1]))
	if not(general_convex_test(l1, l2, p)):
		return False
	return True


def convex_test(l1, final_point, p):
	const_val = 300
	final_points = np.linspace(l1[-1]+0.0001, final_point, const_val)
	ds = (final_point-l1[-1])/const_val # Spacing for each test point.
	for i in range(1, len(final_points)-1):
		diam_mid = diam_p(l1+[final_points[i]], p)
		diam_front = diam_p(l1+[final_points[i-1]], p)
		diam_back = diam_p(l1+[final_points[i+1]], p)
		print((diam_front, diam_back, diam_mid))
		val = diam_front + diam_back - 2*diam_mid
		print(val/(ds**2))
		if val >= 0.000001:
			return False
	return True

def hypothesis_test(l1, new_point, p=3):
	alpha = 1/p
	old_length = l1[-1] - l1[0]
	new_length = new_point - l1[0]
	l2 = l1 + [new_point]
	circum_2_p = diam_p(l2, p)
	circum_1_p = diam_p(l1, p)
	circum_val = circum_2_p-circum_1_p
	length_val = new_length - old_length
	print(circum_val)
	print(length_val)
	print (circum_val - length_val)
	print (circum_val / length_val)
	print ((circum_val - length_val)/length_val)
	return circum_val > length_val - eps

def rad_alpha(ls, alpha=1/3):
	return rad_p(ls, alpha, 1/alpha)

def diam_p(ls, p=3):
	ls = np.sort(ls)
	alpha = 1/p
	return diam_alpha_p(ls, alpha, p)

def diam_alpha_p(ls, alpha=1/3, p = 3):
	return pow(4*rad_squared(ls, alpha), p/2.0)

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
def screw(ls, power):
	n = len(ls)
	M = np.zeros((n, n))
	for i in range(n):
		for j in range(n):
			D = abs(ls[i]-ls[j])
			M[i][j] = pow(D, power)
	return M
