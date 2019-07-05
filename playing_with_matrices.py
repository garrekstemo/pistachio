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
				   
				   
				   
W = np.linalg.eig(Omega)  # eigenvectors
eigval = np.linalg.eigvals(Omega)  # eigenvalues

print(eigval.shape)
eigval.shape = (4,1)
print(eigval.shape)

c = np.dot(W, eigval)

#print(c, "\n")

print(Omega, "\n")
print(eigval)