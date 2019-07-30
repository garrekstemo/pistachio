import numpy as np
import scipy as sp
import scipy.constants as sc
import cmath
import matplotlib.pyplot as plt

c = sc.c

lmbda = 10
n0 = 1.0
n1 = 2.89
n2 = 3.38
a = 0.5*lmbda
b = 0.5*lmbda

omega = 2*np.pi *c / lmbda
theta = [t for t in np.linspace(0.001, np.pi/2-0.001, num=100000)]

def beta_func(theta_, n0_):
	return n0_*omega/c * np.sin(theta_)

def Q(beta_, n_):
	g = (n_*omega / c)
	q_ = cmath.sqrt(beta_**2 - g**2)
	return q_
	
def P(beta_, n_):
	g = (n_*omega / c)
	p_ = cmath.sqrt(g**2 - beta_**2)
	return p_
	
def klx_func(beta_, n_):
	g = (n_*omega / c)
	klx_ = cmath.sqrt(g**2 - beta_**2)
	return klx_


def A_M(q_, p_, a_, b_):
	A = cmath.exp(q_*a_)*( cmath.cos(p_*b_) - 0.5*(p_/q_ - q_/p_)*cmath.sin(p_*b_) )
	return A
	
def B_M(q_, p_, a_, b_):
	B = cmath.exp(-q_*a_)*( -0.5*(p_/q_ + q_/p_)*cmath.sin(p_*b_) )	
	return B
	
def C_M(q_, p_, a_, b_):
	C = cmath.exp(q_*a_)*( 0.5*(p_/q_ + q_/p_)*cmath.sin(p_*b_) )
	return C
	
def D_M(q_, p_, a_, b_):
	D = cmath.exp(-q_*a_)*( cmath.cos(p_*b_) + 0.5*(p_/q_ - q_/p_)*cmath.sin(p_*b_) )
	return D
	
def mode_condition(q_, p_, h_, d):
	cond = 0.5*(1+p_/q_)*cmath.cos(h_*d) - 0.5*cmath.sin(h_*d)*(h_/q_ - p_/h_)
	return cond
	
def chebyshev(A_, D_, N_=1):
	KL = cmath.acos(0.5 * (A_ + D_))
	cheb = A_* cmath.sin(N*KL)/cmath.sin(KL) - cmath.sin((N-1)*KL)/cmath.sin(KL)
	if cheb.imag < 10**-10:
		cheb = cheb.real
	return cheb
	
def construct_matrix(A_, B_, C_, D_, N_=1):
	TM = np.array([[A_, B_], [C_, D_]])
	return TM
	
def E_in(B0_, q0_, x_):
	"""x < 0. Left side"""
	E = B0_ * cmath.exp(q0_*x_)  # Needs to go to zero as x -> -infinity
	return E.real
	
def E_mid(Al, Bl, beta_, klx_, x_, xl):
	"""Valid for x_(l-1) < x < x_l"""
	E = cmath.exp(-1j*klx_*(x_ - xl)) + cmath.exp(1j*klx_*(x_ - xl))
	return E.real

def E_out(qs_, x_, xN):
	"""x > 0. Right side"""
	E = cmath.exp(-qs_*(x_ - xN))
	return E.real

	
modes = []
q0s = []
k1xs = []
k2xs = []
A1 = []
B1 = []

N=1
E0 = np.array([0, 0])
B0 = []
Es = np.array([1, 0])
for t in theta:
	Beta = beta_func(t, n1)
	q0 = Q(Beta, n0)   #beginning
	q = Q(Beta, n1)
	p = P(Beta, n2)
# 	print("Wavenumbers")
# 	print("beta =", Beta)
# 	print("q0 =", q0, n0)  # Must be positive
# 	print("q = ik1x =", q, "k1x =", -1j*q, n1)
# 	print("p = k2x = ", p, n2)

	A = A_M(q, p, a, b)
	B = B_M(q, p, a, b)
	C = C_M(q, p, a, b)
	D = D_M(q, p, a, b)
	cin = chebyshev(A, D)

	if np.round(cin, 4) == 0:
		print("Set for theta = {}".format(np.round(t*180/np.pi, 2)))
		print("    Chebyshev =", cin)
		print("    Beta = {}".format(Beta))
		TM = construct_matrix(A, B, C, D, N)
		E0 = np.matmul(TM, Es)
		B0.append(E0[1])
		TMinv = np.linalg.inv(TM)
		EN = np.matmul(TMinv, E0)
		AN = EN[0]
		BN = EN[1]
		Beta = np.round(Beta, 3).real
		if Beta not in modes:
			modes.append(Beta)
			q0s.append(q0)
			k1xs.append(-1j*q)
			k2xs.append(p)
			A1.append(AN)
			B1.append(BN)

print("\nbeta for guided waves")
print(modes)

KL = int(b+a)

xl = [i for i in np.linspace(-KL, 0, num=1000)]
xm1 = [i for i in np.linspace(0, b, num=1000)]
xm2 = [i for i in np.linspace(b, KL, num=1000)]
xr = [i for i in np.linspace(KL, 2*KL, num=1000)]

E0 = []
E1 = []
E2 = []
Es = []

idx = 0
m = modes[idx]
q0 = q0s[idx]
qs = q0s[idx]
k1x = k1xs[idx]
k2x = k2xs[idx]
print("Testing wavenumbers")
print(p, q)
print("B0 =", B0[idx])

# This currently gives some continuity, but no evanescent waves. UUGGGGH!!
for z in xl:
	E0.append(E_in(B0[idx], q0, z))
for z in xm1:
	E1.append(E_mid(m, k2x, z, b))
for z in xm2:
	E2.append(E_mid(m, k1x, z, a))
for z in xr:
	Es.append(E_out(qs, z, KL))


plt.plot(xl, E0)
plt.plot(xm1, E1)
plt.plot(xm2, E2)
plt.plot(xr, Es)
plt.show()


# Single slab antisymmetric wave guides
# 
# xl = [i for i in np.linspace(-3*d, -d, num=1000)]
# xm  = [i for i in np.linspace(-d, 0, num=1000)]
# xr = [i for i in np.linspace(0, 3*d, num=1000)]
# E3_l = []
# E2_l = []
# E1_l = []
# m = modes[1]
# 
# for i in xr:
# 	E3_l.append(E_field3(m, i))
# 
# for i in xm:
# 	E2_l.append(E_field2(m, i))
# 
# for i in xl:
# 	E1_l.append(E_field1(m, i))
# 	
# plt.plot(xr, E1_l)
# plt.plot(xm, E2_l)
# plt.plot(xl, E3_l)
# plt.show()





