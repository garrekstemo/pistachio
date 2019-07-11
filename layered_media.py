#! /anaconda3/bin python

import csv
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc

c = sc.c  # speed of light
h = sc.h  # planck's constant

class Light:

	def __init__(self, wavelength_):
		
		self.wavelength_ = wavelength_
		self.omega = 2*np.pi*c / wavelength_  # angular frequency
		self.freq = c / wavelength_  # frequency
		self.k = 2*np.pi / wavelength_  # wavenumber
		self.energy_J = h*c / wavelength_  # energy in Joules
		self.energy_ev = self.energy_J / sc.eV  # energy in electronvolts
	

def wavenumber(n, omega, theta=0):
	
	k_x = n*omega/sc.c * np.cos(theta)
	k_z = n*omega/c * np.sin(theta)
	return k_x, k_z


def propagation_matrix(k_, d_):
	"""Inputs: wave number, thickness of medium
	   Output: propagation matrix (phase accumulation for plane wave 
	   			propagating through homogeneous medium)."""	
	phi = k_*d_
	i = 1j
	P_i = np.matrix([[np.exp(i*phi), 0], [0, np.exp(-i*phi)]])
	
	return P_i


def dynamical_matrix(n, theta=0):
	"""Inputs: index of refraction, angle of incidence
		Outputs: dynamical matrices for s-wave and p-wave."""
		
	# s-wave dynamical matrix
	m = n * np.cos(theta)
	Ds = np.matrix([[1, 1], [m, -m]])
# 	Dsinv = np.linalg.inv(Ds)
					 
	# s-wave dynamical matrix
	Dp = np.matrix([[np.cos(theta), np.cos(theta)], [n, -n]])
# 	Dpinv = np.linalg.inv(Dp)
		 
	return Ds, Dp


def make_layer(n, k_x, d, theta=0):

	D = dynamical_matrix(n, theta)[0]
	Dinv = np.linalg.inv(D)
	P = propagation_matrix(k_x, d)
	layer = np.linalg.multi_dot([D, P, Dinv])
	
	return layer
	

def fresnel(n1, n2, k1x, k2x):
	"""Inputs: angular frequency of incident light
			   refractive indices for two media
			   angle of incidence (to calculate wavenumbers
	   Outputs: reflection and transmission coefficients for s and p waves."""
	
	# s-wave reflection, transmission coefficients
	r12_s = (k1x - k2x) / (k1x + k2x)
	t12_s = 2*k1x / (k1x + k2x)
	
	# p-wave reflection, transmission coefficients
	r12_p = (n1*n1*k2x - n2*n2*k1x) / (n1*n1*k2x + n2*n2*k1x)
	t12_p = (2*n1*n1*k2x) / (n1*n1*k2x + n2*n2*k1x)
	
	return r12_s, t12_s, r12_p, t12_p


def transmission_matrix(r_ij, t_ij):
	"""I don't think I'm actually using this yet. Delete comment when used."""
	"""Inputs: Fresnel transmission, reflection coefficients from medium i to medium j.
	   Output: Corresponding transmission matrix linking amplitudes
	   from medium i to medium j."""
	D_ij = 1/t_ij * np.matrix([[1, r_ij], [r_ij, 1]])    
	return D_ij
	

def multilayer_matrix(n0, n1, ns, d, omega, theta=0, num_layers=0):

	TOL = 1e-06  #small number tolerance

	k1x = wavenumber(n1, omega, theta)[0]


	D0 = dynamical_matrix(n0, theta)[0] # Incident medium
	D0inv = np.linalg.inv(D0)
	D1 = dynamical_matrix(n1, theta)[0]
	D1inv = np.linalg.inv(D1)
	Ds = dynamical_matrix(ns, theta)[0] # Transmit medium
	
	L = make_layer(n1, k1x, d)
	
	M = np.linalg.multi_dot([D0inv, L, Ds])

	#M.real[abs(M.real) < TOL] = 0.0
	#M.imag[abs(M.imag) < TOL] = 0.0
	
	# Quarter-wave stack test
# 	M = D0inv*P0*D0 * D1*P1*D1inv  # Quarter-wave stack	
	#M_quarter = np.matrix([[-n1/n0, 0], [0, -n0/n1]])
	
	return M

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


def reflectance(M_):
	"""Input: multilayer matrix, M.
	   Output: reflectance calculation."""
	M21 = M_.item((1, 0))
	M11 = M_.item((0, 0))
	r = M21 / M11
	r_sq = r * np.conj(r)
	r_sq = r_sq.real
	R = r_sq
	   
	return R, r
	   
def transmittance(M_, n0_, ns_, theta_0=0, theta_s=0):
	"""Inputs: Multilayer matrix of dynamical matrices and propagation matrix."""
	M11 = M_.item((0, 0))
	t = 1/M11
	t_sq = t * np.conj(t)
	t_sq = t_sq.real
	
	#T = ((ns_*np.cos(theta_s)) / (n0_*np.cos(theta_0))) * t_sq
	T = np.linalg.det(M_) * t_sq
	
	return T, t
	
	
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

def make_wavelengths(num_points, min_l, max_l):
	"""Make wavelength data from desired number of points, 
		starting and ending wavelengths."""
	wavelen_list = []
	wavelen = min_l
	data_width = (max_l - min_l) / num_points
	
	i = 0	
	while i < num_points:
		wavelen_list.append(wavelen)
		wavelen += data_width
		i+=1
	wavelen_list.sort()
		 
	return wavelen_list
	

def get_wave_data(data_file):
	"""Not Implemented"""
	with open(data_file, 'r') as f:
		wavelength = []
		y = []
		K = []
		reader = csv.reader(f)
		next(reader, None)
		for row in reader:
			wavelength.append(float(row[0]))
			y.append(float(row[1]))
			try:
				K.append(float(row[2]))
			except IndexError:
				pass
	f.close()

	return wavelength, y, K

def main():
	
	# Inputs
	#NUM_PTS = 100
	#MIN = 200  #nm
	#MAX = 20000  #nm
	n_air = 1.0
	n_other = 1.4
	d = 35*10**-9  # thickness of medium
	inc_ang = 0.

	data_Au_n = '/Users/garrek/projects/pistachio/data/ag_n.csv'
	data_Au_T = '/Users/garrek/projects/pistachio/data/transmittance_au_SiO2_Ciesielski.csv'

	# Get a bunch of downloaded data
	wavelen_n, Au_n, Au_K = get_wave_data(data_Au_n)
	wavelen_T = get_wave_data(data_Au_T)[0]  # experimental data for Au T
	Au_T = get_wave_data(data_Au_T)[1]
	
	# Outputs
	R = []  # my reflectance data
	T = []  # my transmittance data
	
	NUM_POINTS = len(Au_n)
	i = 0
	while i < NUM_POINTS:
# 		print("="*60)
		lmbda = wavelen_n[i] * 10**-6
		n_au = Au_n[i] - 1j*Au_K[i]
		light = Light(lmbda)
		omega = light.omega

		#Get wavenumbers for two media
		k0x = wavenumber(n_air, omega)[0]
		k1x = wavenumber(n_au, omega)[0]
		
		#Use dynamical matrix
		M_i = multilayer_matrix(n_air, n_au, n_other, d, omega)  
				
		#Use Fresnel equations
		r_01 = fresnel(n_air, n_au, k0x, k1x)[0]
		t_01 = fresnel(n_air, n_au, k0x, k1x)[1]
		D_01 = transmission_matrix(r_01, t_01)
		r_10 = fresnel(n_au, n_air, k1x, k0x)[0]
		t_10 = fresnel(n_au, n_air, k1x, k0x)[1]
		D_10 = transmission_matrix(r_10, t_10)

		trans = transmittance(M_i, n_air, n_other)[0].real
		refl = reflectance(M_i)[0]
		T.append(trans)
		R.append(refl)

		i+=1
	
	
	#Make Plots
	fig, axs = plt.subplots(3, 1, sharex=True)
	
	ax = axs[0]
	ax.plot(wavelen_n, Au_n, label="n, downloaded data")
	ax.plot(wavelen_n, Au_K, label="K, downloaded data")
	ax.set_ylabel('refractive index for Au')
	ax.legend()
	
	ax = axs[1]
	ax.plot(wavelen_T, T, label="calculated data")
	ax.plot(wavelen_T, Au_T, label="downloaded data")
	ax.set_ylabel('transmittance')
# 	ax.set_xlim(0, 5)
	ax.legend()
	
	ax = axs[2]
	ax.plot(wavelen_T, R, label="calculated data")
	ax.set_xlabel('wavelength ($\mu$m)')
	ax.set_ylabel('reflectance')
	ax.legend()
	
	fig.tight_layout()
	plt.show()
	

if __name__ == '__main__':
	main()
