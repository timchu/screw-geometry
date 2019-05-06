import screw as sc
import numpy as np

# A list of numbers a_1, ... a_n, followed by np1 < np2.
""" Test whether 
R(a_1. ... a_n, np2)^1/alpha - R(a_1, ... a_n, np1^1/alpha)
is bigger than the same thing without a_n
"""
def incremental_inequality(ls, np1, np2, alpha=1/10):
	l1 = ls + [np1]
	l2 = ls + [np2]
	ls1 = ls[:-1] + [np1]
	ls2 = ls[:-1] + [np2]
	V1 = sc.rad_alpha(l2, alpha) - sc.rad_alpha(l1, alpha)
	V2 = sc.rad_alpha(ls2, alpha) - sc.rad_alpha(ls1, alpha)
	length1_alpha = pow(pow(np1-ls[0], alpha)/2, 1/alpha)
	length2_alpha = pow(pow(np2-ls[0], alpha)/2, 1/alpha)
	V3 = length2_alpha - length1_alpha
	print(V1)
	print(V2)
	print(l2)
	print(l1)
	print(ls1)
	print(ls2)
	print(V1/V2)
	print(V1/V3)
	return V1 / V2 > 1-sc.eps

def other_elem_removal_ii(ls, np1, np2, alpha=1/10, removal=2):
	l1 = ls + [np1]
	l2 = ls + [np2]
	ls1 = ls[:-removal] + [np1]
	ls2 = ls[:-removal] + [np2]
	V1 = sc.rad_alpha(l2, alpha) - sc.rad_alpha(l1, alpha)
	V2 = sc.rad_alpha(ls2, alpha) - sc.rad_alpha(ls1, alpha)
	print(V1)
	print(V2)
	print(V1/V2)
	return V1 / V2 > 1-sc.eps
