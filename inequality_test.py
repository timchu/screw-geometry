import numpy as np
from numpy.linalg import inv
from screw import screw
""" Tested a number of inequalities. Starting with:
 Whether for Euclidean distance matrices (hereafter referred to as distance matrices):
 1. Radius(D1)+Radius(D2) <= 2*Radius(D1 + D2 / 2). 
 	Results: False (dist_inequality)
 2. Whether a^T A a + b^T B b <= 2*(a+b)/2 ^T (A+B)/2 (a+b)/2
 	Results: False (big_concavity_test)
 3. Tested some p-norm inequality for diagonal matrices. See val, inequality for specs.
 4. Tested whether distance matrices were preserved under harmonic mean. This was false. (See Harmonic, neg_type)
 """

def rand_3_dist():
	return n_rand_dist(3)
def n_rand_dist(n):
	# Tends to generate things like a simplex.
	a = np.random.rand(n,n-1)
	A = np.zeros((n,n))
	for i in range(n):
		for j in range(n):
			A[i][j] = np.linalg.norm(a[i]-a[j])**2
	return A
def radius_thing(D):
	n = len(D)
	ones = np.ones(n)
	two_rad_squared = 1/(ones @ np.linalg.inv(D) @ ones) # You need to invert to get the radius
	b = np.linalg.inv(D)@ones*two_rad_squared
	print ("Barycentric Coords are:")
	return two_rad_squared

# concavity Test:
# Test whether a^T A a + b^T B b <= 2* (a+b)/2^T (A+B/2)(a+b)/2
# Result: fails, for A and B from distance matrix |i-j|. See big_concavity_test
# Even the hacked concavity test fails.
def hacked_concavity_test(a,A,b,B):
	concavity_vals = [(c*a + (1-c)*b) @ (A+B)/2 @ (c*a+(1-c)*b) for c in np.linspace(-10, 10, 3000)]
	return  a @ A @ a + b @ B @ b - 2 * max(concavity_vals)

def concavity_test(a, A, b, B):
	concavity_val =  a @ A @ a + b @ B @ b - 2 * (a+b)/2 @ (A+B)/2 @ (a+b)/2
	return concavity_val

# Even for screw distance with edge-squared metric, Fails this concavity test.
def big_concavity_test():
	for i in range(100):	
		ls = np.sort(np.random.rand(3))
		A = screw(ls, 1)
		ls2 = np.sort(np.random.rand(3))
		B = screw(ls2, 1)
		(a,b) = np.random.rand(3), np.random.rand(3)
		(a,b) = a/sum(a), b/sum(b)
		print(i)
		if (concavity_test(a, A, b, B) >= 0.00001):
			print("A: ", A)
			print("a: ", a)
			print("B: ", B)
			print("b: ", b)
			print("Avg: ", (A+B)/2)
			print("a+b/2: ", (a+b)/2)
			print(a @ A @ a)
			print(b @ B @ b)
			print((a+b)/2 @ (A+B)/2 @ (a+b)/2)
			print("Concavity Value:", concavity_test(a,A,b,B))
			print("Fails second concavity test")
			return False
	return True
# New Concavity Test:
# Test whether a^T A a + b^T B b <= avg thing
# when A and B are the same matrix but with only the last distance changed.
def new_big_concavity_test():
	for i in range(100):	
		n = 3
		ls = list(np.sort(np.random.rand(n)))
		A = screw(ls, 1)
		ls2 = ls[:-1] + [ls[-1]+np.random.rand()]
		B = screw(ls2, 1)
		(a,b) = np.random.rand(n), np.random.rand(n)
		(a,b) = a/sum(a), b/sum(b)
		print(i)
		if (concavity_test(a, A, b, B) >= 0.00001):
			print(ls, ls2, a,b)
			print("")
			print(A,B,(A+B)/2)
			print("Concavity Value:", concavity_test(a,A,b,B))
			print("Fails big concavity test")
			return False
	return True

def p_hacked_concavity_test():
	for i in range(100):	
		n = 3
		ls = list(np.sort(np.random.rand(n)))
		A = screw(ls, 1)
		ls2 = ls[:-1] + [ls[-1]+np.random.rand()]
		B = screw(ls2, 1)
		(a,b) = np.random.rand(n), np.random.rand(n)
		(a,b) = a/sum(a), b/sum(b)
		print(i)
		if (hacked_concavity_test(a, A, b, B) >= 0.00001):
			print(ls, ls2, a,b)
			print("")
			print(A,B,(A+B)/2)
			print("Concavity Value:", hacked_concavity_test(a,A,b,B))
			print("Fails hacked concavity test")
			return False
	return True

# I attempted to test whether if D is negative type, then 
# max(x, y, z, orthogonal to ones) x^T D1 x + y^T D2 y <= z^T avg(D1, D2) z
# It seems always false :D :D :D :D :D :D :D :D.
def dist_inequality(p):
	D1 = rand_3_dist()
	D2 = rand_3_dist()
	D3 = (D1+D2)/2
	print(D1)
	print(D2)
	print(D3)
	rad_D1 = radius_thing(D1)
	rad_D2 = radius_thing(D2)
	rad_D3 = radius_thing(D3)
	print (rad_D1+rad_D2-2*rad_D3)
	return (rad_D1+rad_D1-2*rad_D3 <= 0.000001)

# Test for diagonal matrix with two entries.
# It seems to work.
def val(a, b, p):
	av = pow(a, 1/p)
	bv = pow(b, 1/p)
	# Minimal value of ax^2 + by^2 when x+y=1
	v = (2*av*bv)/(av+bv)
	# Inverted powered v should be the circumradius
	return pow(1/v, p)

def inequality(a, b, c, d, p = 3/2):
	firstval = val(a,b,p)
	secondval = val(c, d, p)
	thirdval = val((a+c)/2, (b+d)/2, p)
#	print(firstval + secondval)
#	print (2*thirdval)
#	print("Values are:")
#	print(firstval, secondval, thirdval)
	return firstval + secondval - 2*thirdval
def ineq_test():
	for _ in range(1000):
		a = np.random.rand()+0.001
		b = np.random.rand()+0.001
		c = np.random.rand()+0.001
		d = np.random.rand()+0.001
		p = 10*np.random.rand()+0.1
		p=1
		s = inequality(a,b,c,d,p)
		print(s)
		if s <  -0.0000001:
			print((a,b,c,d,p))
			return False
	print(p)
	return True

# Test for 3-entry diagonal matrix.
def val_3(a, b, c, p):
	av = pow(a, 1/p)
	bv = pow(b, 1/p)
	cv = pow(c, 1/p)
	# Minimal value of ax^2 + by^2 + cz^2 when x+y+z = 1 
	# ax = by = cz, x = 1/a, y = 1/b, z = 1/c, normalized
	# value =1/(1/a+1/b+1/c)^2*(1/a+1/b+1/c)
	v = 1/(1/a+1/b+1/c)
	return pow(1/v, p)
def inequality_3(a,b,c,d,e,f,p=3/2):
	firstval = val_3(a,b,c,p)
	secondval = val_3(d,e,f, p)
	thirdval = val_3((a+d)/2, (b+e)/2, (c+f)/2, p)
#	print(firstval + secondval)
#	print (2*thirdval)
#	print("Values are:")
#	print(firstval, secondval, thirdval)
	return firstval + secondval - 2*thirdval

def three_ineq_test():
	for _ in range(1000):
		a = np.random.rand()+0.001
		b = np.random.rand()+0.001
		c = np.random.rand()+0.001
		d = np.random.rand()+0.001
		e = np.random.rand()+0.001
		f = np.random.rand()+0.001
		p = 4*np.random.rand()+0.5
		q = 1
		s = inequality_3(a,b,c,d,e,f,p)
		sq = inequality_3(a,b,c,d,e,f,q)
		t = inequality_3(a,a,a,a,a,b,p)
		tq = inequality_3(a,a,a,a,a,b,q)
		u = inequality_3(a,b,a,a,c,b,p)
		uq = inequality_3(a,b,a,a,c,b,q)
		if s <  -0.0000001:
			print((a,b,c,d,e,f,p))
			print(s)
			return False
		if sq <  -0.0000001:
			print((a,b,c,d,e,f,q))
			print(sq)
			return False
		if t < 0.0001:
			print((a,a,a,a,a,b,p))
			print("Failed on t")
			print(t)
			return False
		if tq < -0.00000001:
			print((a,a,a,a,a,b,q))
			print("Failed on tq")
			print(tq)
			return False
		if u < -0.00000001:
			print((a,a,a,a,a,b,p))
			print("Failed on u")
			return False
		if uq < -0.00000001:
			print((a,a,a,a,a,b,q))
			print("Failed on uq")
			return False
	print(p)
	print(q)
	return True

# Testing whether negative type matrices are preserved under harmonic mean.	
# They probably aren't. If the experiments show it, I'm probably not being thorough enough.

def harmonic_mean(D1, D2):
	return 2*inv(inv(D1)+inv(D2))
def neg_type(D):
	n = len(D)
	P = np.zeros((n-1, n-1))
	for i in range(1, n):
		for j in range(1, n):
			P[i-1][j-1] = D[0][i]+D[j][0]-D[i][j]
	print( np.linalg.eig(P)[0])
	sorted_real_eigvals =  np.sort(np.real(np.linalg.eig(P)[0]))
	print(sorted_real_eigvals[0])
	return sorted_real_eigvals[0] >-.000001


def structured_n_dist(n):
	x_size = int(n/2)
	y_size = n-x_size
	dim = n
	x = np.random.rand(x_size,dim)
	y = np.random.rand(y_size,dim)
	for i in range(len(y)):
		y[i][0] += 2
	a = np.concatenate((x,y))
	A = np.zeros((n,n))
	for i in range(n):
		for j in range(n):
			A[i][j] = np.linalg.norm(a[i]-a[j])**2
	return A
