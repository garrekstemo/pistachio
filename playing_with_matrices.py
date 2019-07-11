#!/anaconda3/bin python

import numpy as np

eps_r = 9.
mu_r = 1.
n = np.sqrt(eps_r * mu_r)

# Normal incidence (no x, y component for wave number)
kx = 0
ky = 0



Omega = np.matrix([[0, 0, (kx*ky / eps_r), (mu_r - kx*kx / eps_r)],
				   [0, 0, (ky*ky / eps_r - mu_r), (-kx*ky / eps_r)],
				   [(kx*ky / mu_r), (eps_r - kx*kx / mu_r), 0, 0],
				   [(ky*ky/mu_r - eps_r), (-kx*ky/mu_r), 0, 0]])
				   
P = 1/eps_r * np.matrix([
			   [kx*ky, (mu_r*eps_r - kx*kx)],
			   [(ky*ky - mu_r*eps_r), -kx*ky]
			   ])
				   
				   
Q = 1/mu_r * np.matrix([
						[kx*ky, (mu_r*eps_r - kx*kx)],
						[(ky*ky - mu_r*eps_r), -kx*ky]
						])
						
# P.Q is the same as Omega.Omega, but with different dimensions, since P and Q operate on 2,1 vectors.					
print("P Matrix:\n", P)
print("Q Matrix:\n", Q)
print("P . Q:\n", np.dot(P,Q))
print("Omega . Omega:\n", np.dot(Omega, Omega))

				   
W, V = np.linalg.eig(Omega)  # linalg.eig produces both eigenvectors and eigenvalues
#V = np.linalg.eigvalsh(Omega)  # eigenvalues (hermitian if there's an h)

c = np.dot(W, V)

print("Eigenvectors")
print(W)
print("Eigenvalues")
print(V)

print("W.V")

print(c)
