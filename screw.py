import numpy as np
from scipy.linalg import fractional_matrix_power as frac_power
eps = pow(10.0, -7)

def pinv_partial_frac_lap(alpha=1/3.0, num_points = 100, partial = 8, inv = True):
	lap = np.zeros((num_points, num_points))
	for i in range(len(lap)):
		lap[i][i] = 2
		if i >= 1:
			lap[i-1][i]=-1
			lap[i][i-1]=-1
	lap[0][0]=1
	lap[-1][-1]=1
	print(lap)
#	lap = lap/(partial)
	if (inv):
		fractional_lap = np.real(frac_power(np.linalg.pinv(lap),0.5+alpha))
	else:
		fractional_lap = np.real(frac_power(lap, 0.5+alpha))
	start = int(num_points/2-partial/2)
	M = fractional_lap[start:start+partial, start:start+partial]
	print(np.around(M, 2))
	return M

def partial_fractional_lap(alpha=1/3.0, num_points = 100, partial=8, inv=True):
	lap = np.zeros((num_points, num_points))
	for i in range(len(lap)):
		lap[i][i] = 2
		if i >= 1:
			lap[i-1][i]=-1
			lap[i][i-1]=-1
	lap[0][0]=2
	lap[-1][-1]=1
#	lap = lap/(partial)
	if (inv):
		fractional_lap = np.real(frac_power(lap, -0.5-alpha))
	else:
		fractional_lap = np.real(frac_power(lap, 0.5+alpha))
	start = int(num_points/2-partial/2)
	M = fractional_lap[start:start+partial, start:start+partial]
	print(np.around(M, 2))
	return M

def inner_prod(alpha=1/3.0, num_points=1000, partial=200, inv=True):
	M = np.zeros((partial+1, partial))
	for x in range(partial):
		M[0][x] = -	1
		M[x+1][x] = 1
	S =  np.dot(np.dot(M.T, partial_fractional_lap(alpha, num_points, partial+1, inv)), M)
	return S
# DEPRECATED: Replace power with alpha before using.
def comparator(i, j, alpha=1/3.0):
	return (pow(i, 2*alpha) + pow(j, 2*alpha) - pow(abs(i-j), 2*alpha))/2

# DEPRECATED: Replace power with alpha before using.
def fake_screw_circumcenter_barycentric(alpha=1/3.0, num_points=1000, interval=1, line_size=10):
	n = int(num_points*interval/line_size)
	lap = np.zeros((num_points, num_points))
	partial_frac_lapl = partial_fractional_lap(alpha, num_points, n, False)
	vec = np.zeros(n)
	for i in range(n):
		vec[i] = pow(i+1, 2*alpha)
#	print(vec)
#	print(lap)
#	print(np.around(np.dot(fractional_lap, fractional_lap)),3)
	#print(fractional_lap)
	#print(partial_frac_lap)
	return np.around(np.dot(partial_frac_lap, vec), 3)

# DEPRECATED: Replace power with alpha before using.
def circumcenter_barycentric_screw(alpha=1/3.0, num_points=4, interval=1):
	n = num_points
	M = np.zeros((n, n))
	for i in range(n):
		for j in range(n):
			M[i][j] = (pow(i+1, 2*alpha) + pow(j+1, 2*alpha) - pow(abs(i-j), 2*alpha))
	vec = np.zeros(n)
	for i in range(n):
		vec[i] = pow(i+1, 2*alpha)
	#print(M)
	#print(vec)
	return np.around(np.linalg.solve(M, vec), 3)

# DEPRECATED: Replace power with alpha before using.
def screw_segment(power=2.0/3, num_points=2200):
	interval_size = 1
	l1 = np.linspace(0, interval_size, num_points)
	q = rad_squared(l1, power)
	print(q)
	return q

# DEPRECATED: Replace power with alpha before using.
def line_test(back=0, numpoints=800, weak=1, gap=1, power=2.0/3):
	interval_size = 10
	l1 = np.linspace(back, interval_size, numpoints)
	new_point = l1[-1]+interval_size*gap
	if (numpoints < 10):
		print(l1)
	print(new_point)
	if (weak == 1):
		if not(weak_hypothesis_test(l1.tolist(), new_point, power)):
			return False
	else:
		if not(hypothesis_test(l1.tolist(), new_point, power)):
			return False
	return True

# DEPRECATED: Replace power with alpha before using.
def distorted_gap_hypothesis_test(numpoints=30, gap=1,
	distortion_factor=0.1, power=2.0/3):
	def inject_noise(ls, distortion_factor):
		l1 = np.copy(ls)
		for i in range(1, len(l1)-1):
			l1[i] += (l1[i]-l1[i-1])*np.random.rand()*distortion_factor
		return l1
	interval_size = 10
	if distortion_factor > -eps:
		l1 = inject_noise(np.linspace(0, interval_size, numpoints),
			distortion_factor)
	elif distortion_factor == -1:
		l1 = interval_size*np.sort(np.random.rand(numpoints))
		l1[0]=0
		l1[-1]=interval_size
	elif distortion_factor == -2:
		l1 = np.linspace(0, 2*interval_size/3, numpoints)
		l1[-1]=interval_size
	elif distortion_factor==-3:
		l1 = np.linspace(interval_size/3,
			interval_size, numpoints)
		l1[0]=0
	print(l1)
	new_point = l1[-1]+interval_size*gap
	print(new_point)
	if not(hypothesis_test(l1.tolist(), new_point, power)):
		return False
	return True
# DEPRECATED: Replace power with alpha before using.
def fixed_gap_hypothesis_test(numpoints=30, gap=3,
	power=2.0/3):
	interval_size = 100
	l1 = np.linspace(0, interval_size, numpoints)
	print(l1)
	new_point = l1[-1]+interval_size*gap
	print(new_point)
	if not(hypothesis_test(l1.tolist(),
		new_point, power)):
		return False
	return True

def mass_hypothesis_test(n=1000, alpha=1.0/3, numpoints=10):
	# Generate a bunch of points on a line to test this.
	for i in range(n):
		test_numpoints = np.random.randint(1, numpoints)
		l1 = np.sort(test_numpoints *
			np.random.rand(test_numpoints))
		new_point = l1[-1]+np.random.rand()
		if not(hypothesis_test(l1.tolist(), new_point, alpha)):
			return False
	return True

# Add new point to the end, see how it effects the Circumradius_2^p-Circumradius_1^p quantity.
# I want to test whether [-1, 0, 100] has a larger p_circum_diff or [0, 100]
def p_circumradius_change(alpha,l):
	l2 = list(np.sort(l))
	l1 = l2[:-1]
	circum_2_p = pow(rad_squared(l2, alpha), 1.0/(2*alpha))
	circum_1_p = pow(rad_squared(l1, alpha), 1.0/(2*alpha))
	print(circum_2_p)
	print(circum_1_p)
	print(circum_2_p-circum_1_p)
	print("")
	return circum_2_p - circum_1_p
def square_circumradius_change(alpha, l):
	l2 = list(np.sort(l))
	l1 = l2[:-1]
	circum_2_p = rad_squared(l2, 2*alpha)
	circum_1_p = rad_squared(l1, 2*alpha)
	print(circum_2_p)
	print(circum_1_p)
	print(circum_2_p-circum_1_p)
	print("")
	return circum_2_p - circum_1_p

def hypothesis_test(l1, new_point, alpha=1.0/3, p = 2):
	old_length = l1[-1] - l1[0]
	new_length = new_point - l1[0]
	l2 = l1 + [new_point]
	circum_2_p = pow(rad_squared(l2, alpha), p/2.0)
	circum_1_p = pow(rad_squared(l1, alpha), p/2.0)
	circum_val = circum_2_p-circum_1_p
	new_val = pow(pow(new_length, alpha)/2, p) - pow(pow(old_length, alpha)/2, p)
	print(circum_val)
	print(new_val)
	print(circum_val/new_val)
	return circum_val > new_val - eps

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
