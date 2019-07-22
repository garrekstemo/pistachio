#! /anaconda3/bin python

import numpy as np
from transfer_matrix import *


def testing_quarterwave(n0=1.0, n1=2.32, n2=1.36, ns=2.32):
	print("="*70)
	print("TESTING AGAINST YEH DATA PG 111-112 -- Result: it works.")
	d = 1
	k1x = np.pi/2
	num_layers=0
	
	D0 = dynamical_matrix(n0)[0] # Incident medium
	D0inv = np.linalg.inv(D0)
	Ds = dynamical_matrix(ns)[0] # Transmit medium
	L1 = make_layer(n1, k1x, d)
	L2 = make_layer(n2, k1x, d)
	
	L = np.dot(L2, L1)
# 	print("This is L: \n", L)
		
	L = np.linalg.matrix_power(L, 1)
# 	print("L^{}".format(num_layers), L)

	# For N=0,1,2,3,4,5
	M0 = np.linalg.multi_dot([D0inv, Ds])
	M1 = np.linalg.multi_dot([D0inv, L1,L2, Ds])
	M2 = np.linalg.multi_dot([D0inv, L1,L2,L1,L2, Ds])
	M3 = np.linalg.multi_dot([D0inv, L1,L2,L1,L2,L1,L2, Ds])
	M4 = np.linalg.multi_dot([D0inv, L1,L2,L1,L2,L1,L2,L1,L2, Ds])
	M5 = np.linalg.multi_dot([D0inv, L1,L2,L1,L2,L1,L2,L1,L2,L1,L2, Ds])
	
	R0 = reflectance(M0)[0]
	R1 = reflectance(M1)[0]
	R2 = reflectance(M2)[0]
	R3 = reflectance(M3)[0]
	R4 = reflectance(M4)[0]
	R5 = reflectance(M5)[0]
	
	R = [R0, R1, R2, R3, R4, R5]
	
	i=0
	while i< len(R):
		print("N", i, ", R =", R[i])
		i+=1
	return 0


def conservation_tests(M_, n0_, ns_, theta_0=0, theta_s=0):
	"""Test conservation laws for lossless media"""
	print("_"*50)
	print("Checking conservation principles\n")
	M11 = M_.item((0, 0))
	M12 = M_.item((0, 1))
	M21 = M_.item((1, 0))
	M22 = M_.item((1, 1))
	MM = np.linalg.det(M_)
# 	print("det M = ", MM)
	MMs = np.conj(MM)
	
	MM = np.sqrt(MM * MMs)
	
	r = M21 / M11
	r_p = - M12 / M11
	t = 1/M11
	t_p = MM / M11
	t_sq = t * np.conj(t)
	t_p_sq = t_p * np.conj(t_p)
	
	R = r * np.conj(r)
	T = ns_ * np.cos(theta_s) / n0_ * np.cos(theta_0) * t_sq
	T_p = n0_ * np.cos(theta_0) / ns_ * np.cos(theta_s) * t_p_sq
	
	# Conditions
	E = R + T
	
	test_1 = t * np.conj(t_p) + r*np.conj(r)  #equals 1
	test_1 = np.round(test_1, 2)
	test_2 = t*np.conj(r_p) + r*np.conj(t)  # equal 0
	test_2 = np.round(test_2, 2)
	test_3 = MM * t  # equal t'
	
	if np.round(E, 1) != 1:
		print("The universe broke.")
		print("R + T = ", R, " + ", T, " = ", E)
	print("\nCheck energy conservation, R+T=1, using t and r.")
	print("T = {}".format(T.real))
	print("R = {}".format(R.real))
	print("E = {}".format(E.real))
	print("\nCheck reflection and transmission coefficients")
	print("tt'* + rr* = ", test_1, "  should be 1")
	print(test_1.real == 1)
	print("tr'* + rt* = ", test_2, "  should be 0")
	print(test_2.real == 0)
	print("\nCheck that T=T'\n", T, " = ", T_p)
	print(np.round(T, 3) == np.round(T_p, 3))
	print("\nCheck that det(M)t = t'")
	print("t' = ", t_p, " = ", test_3)
	print(np.round(t_p, 3) == np.round(test_3, 3))
	print("\nDone with conservation")
	print("_"*50)
