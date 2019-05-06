import numpy as np

# Returns a rectangle with elements +-ls
def rect(ls):
	n = len(ls)
	if (n == 1):
		only_elem = ls[0]
		return np.array([[only_elem, -only_elem]])
	old_rect = rect(ls[:-1])
	num_old_points = int(2**(n-1))
	bot = ls[-1]*np.ones((1, num_old_points))
	M_1 = np.concatenate((old_rect, bot), axis=0)
	M_2 = np.concatenate((old_rect,-bot), axis=0)
	M = np.concatenate((M_1, M_2), axis=1)
	print(M)
	return M
# Columns of P are points
def dist_squared(P):
	n = len(P[0])
	D = np.zeros((n, n))
	for i in range(n):
		for j in range(n):
			print(P[:, i])
			print(P[:, j])
			D[i][j] = np.linalg.norm(P[:,i]-P[:, j])**2
	return D

def rect_dist_squared(ls):
	return dist_squared(rect(ls))






