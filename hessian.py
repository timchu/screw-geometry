import screw as s
import numpy as np

eps = 0.001

def diam_p(ls):
	return s.diam_p(ls, 4)

def mklist(gaps):
	n = len(gaps)+1
	ls = [0 ] * n
	for i in range(1, n):
		ls[i] = gaps[i-1]+ls[i-1]
	return ls

def approx_hessian(func, ls):
	# We're finding the hessian w.r.t gaps, and there are n-1 gaps.
	n = len(ls)
	A = np.zeros((n-1,n-1))
	gaps = [ls[i]-ls[i-1] for i in range(1, n)]
	for i in range(n-1):
		for j in range(n-1):
			ls_1 = ls.copy()
			ls_1[i] += eps
			ls_2 = ls.copy()
			ls_2[j] += eps
			ls_3 = ls.copy()
			ls_3[i] += eps
			ls_3[j] += eps
			A[i][j] = ((func(ls_3)- func(ls_2))/eps-(func(ls_1) - func(ls))/ eps)/eps
	print(ls)
	print(A)
	return A
def gaps_approx_hessian(func, ls):
	# We're finding the hessian w.r.t gaps, and there are n-1 gaps.
	n = len(ls)
	A = np.zeros((n-1,n-1))
	gaps = [ls[i]-ls[i-1] for i in range(1, n)]
	for i in range(n-1):
		for j in range(n-1):
			gaps_1 = gaps.copy()
			gaps_1[i] += eps
			gaps_2 = gaps.copy()
			gaps_2[j] += eps
			gaps_3 = gaps.copy()
			gaps_3[i] += eps
			gaps_3[j] += eps
			ls_1 = mklist(gaps_1)
			ls_2 = mklist(gaps_2)
			ls_3 = mklist(gaps_3)
			A[i][j] = ((func(ls_3)- func(ls_2))/eps-(func(ls_1) - func(ls))/ eps)/eps
	print(A)
	return A
	# f(x+i+j)-f(x+i)-f(x+i)+f(x)
	# Deriv = lim(lim f(x+i+j)-f(x+j)/i) - f(x+j)-f(x)

#INconclusive test. For some epsilon, the test doesn't work when p=2.

def convexity_test():
	for i in range(1000):
		n  = 3 
		# Ensures the average spacing is around 10.
		ls = 10*n*np.sort(np.random.rand(n))
		# Ensures that two points aren't closer together than epsilon.
		ls = [ls[i] + 3*i*eps for i in range(ls)]
		if not(NSD(approx_hessian(diam_p, ls))):
				return False
	return True


def NSD(P):
	return(PSD(-P))

def PSD(P):
	n = len(P)
#	print( np.linalg.eig(P)[0])
#	print((np.linalg.eig(P)[1].T)[0])
#	print(np.linalg.eig(P)[1])
	sorted_real_eigvals =  np.sort(np.real(np.linalg.eig(P)[0]))
#	print(sorted_real_eigvals[0])
# The second condition is some condition number phenom. The condition number is huge for most tests.
	if sorted_real_eigvals[0] <-.00001:# and -sorted_real_eigvals[-1]/sorted_real_eigvals[0] < 100000:
		print(sorted_real_eigvals[0])
		print(-sorted_real_eigvals[-1]/sorted_real_eigvals[0])
		return False
	return True
	# Note: This breaeks for me hwen my threshold is 0.000001


	

