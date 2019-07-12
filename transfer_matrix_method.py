#! /anaconda3/bin python

# Convention used
# Psi(x, t) = Psi_0 * exp(i(kx - wt))
# n = n' + i*n'', where n' is real part of refractive index and n'' is imaginary part.

import csv
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import TMM_tests as TMM

c = sc.c  # speed of light
h = sc.h  # planck's constant


class Light:

	def __init__(self, wavelength):
		
		self.wavelength_ = wavelength
		self.omega = 2*np.pi*c / wavelength  # angular frequency
		self.freq = c / wavelength  # frequency
		self.k = 2*np.pi / wavelength  # wavenumber
		self.energy_J = h*c / wavelength  # energy in Joules
		self.energy_ev = self.energy_J / sc.eV  # energy in electronvolts


class Layer:

	def __init__(self, material, index, thickness):
		self.material = material
		self.index = index
		self.thickness = thickness

	def wavenumber(self, omega, theta=0):
		"""Outputs the wavenumbers for the dielectric for the given
		   angular frequency and angle"""
		k_x = self.index*omega/sc.c * np.cos(theta)
		k_z = self.index*omega/sc.c * np.sin(theta)
		return k_x, k_z

	def propagation_matrix(self, wavenumber):
		"""Inputs: wave number, thickness of medium
		   Output: propagation matrix (phase accumulation for plane wave 
					propagating through homogeneous medium)."""	
		phi = wavenumber*self.thickness
		P_i = np.matrix([[np.exp(-1j*phi), 0], [0, np.exp(1j*phi)]])
	
		return P_i

	def dynamical_matrix(self, theta=0.):
		"""Inputs: index of refraction, angle of incidence
			Outputs: dynamical matrices for s-wave and p-wave."""
		n = self.index
		
		# s-wave dynamical matrix
		m = n * np.cos(theta)
		Ds = np.matrix([[1, 1], [m, -m]])
					 
		# s-wave dynamical matrix
		Dp = np.matrix([[np.cos(theta), np.cos(theta)], [n, -n]])
		
		return Ds, Dp

	def structure(n, k_x, d, theta=0):
		"""Maybe deprecate this"""
		k_x = self.wavenumber()[0]
		D = dynamical_matrix(self.index, self.theta)[0]
		Dinv = np.linalg.inv(D)
		P = propagation_matrix(k_x, self.thickness)
		layer = np.linalg.multi_dot([D, P, Dinv])
	
		return layer
		
		
def multilayer_matrix(array):
	TM = np.linalg.multi_dot(array)
	return TM

def inverse(M):
	Minv = np.linalg.inv(M)
	return Minv

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
	"""Inputs: Fresnel transmission, reflection coefficients from medium i to medium j.
	   Output: Corresponding transmission matrix linking amplitudes
	   from medium i to medium j."""
	D_ij = 1/t_ij * np.matrix([[1, r_ij], [r_ij, 1]])    
	return D_ij


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

	data_Au_n = '/Users/garrek/projects/pistachio/data/Au_SiO2_n_Ciesielski.csv'
	data_Au_T = '/Users/garrek/projects/pistachio/data/transmittance_au_SiO2_Ciesielski.csv'

	# Get a bunch of downloaded data
	wavelen_n, Au_n, Au_K = get_wave_data(data_Au_n)
	wavelen_T = get_wave_data(data_Au_T)[0]  # experimental data for Au transmission
	Au_T = get_wave_data(data_Au_T)[1]
	
	# Outputs
	R = []  # my reflectance data
	T = []  # my transmittance data
	
	NUM_POINTS = len(Au_n)
	i = 0
	while i < NUM_POINTS:

		lmbda = wavelen_n[i] * 10**-6
		n = Au_n[i] + 1j*Au_K[i]
		light = Light(lmbda)
		omega = light.omega
		air = Layer('air', 1.0, 0.)
		au = Layer('gold', n, d)
		other = Layer('other', n_other, 0.)

		#Get wavenumbers for two media
		k0x = air.wavenumber(omega)[0]
		k1x = au.wavenumber(omega)[0]
		
		#Use dynamical matrix
		D0 = air.dynamical_matrix()[0]
		D0inv = inverse(D0)
		D1 = au.dynamical_matrix()[0]
		D1inv = inverse(D1)
		P1 = au.propagation_matrix(k1x)
		D2 = other.dynamical_matrix()[0]
		
		M_list = [D0inv, D1, P1, D1inv, D2]
		M_i = multilayer_matrix(M_list)

		trans = transmittance(M_i, air.index, other.index)[0].real
		refl = reflectance(M_i)[0]
		T.append(trans)
		R.append(refl)

		i+=1
	
# 	TEST = TMM.conservation_tests(M_i, n_air, n_other)
# 	Quarter = TMM.testing_quarterwave()
	
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
# 	plt.show()
	

if __name__ == '__main__':
	main()
