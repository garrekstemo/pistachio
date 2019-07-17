#! /anaconda3/bin python

# Convention used
# Psi(x, t) = Psi_0 * exp(i(kx - wt))
# n = n' + i*n'', where n' is real part of refractive index and n'' is imaginary part.

import argparse
import csv
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.constants as sc
import scipy.interpolate
import TMM_tests as TMM
from ruamel_yaml import YAML

c = sc.c  # speed of light
h = sc.h  # planck's constant
yaml = YAML()


class Light:

	def __init__(self, wavelength):
		
		self.wavelength_ = wavelength
		self.omega = 2*np.pi*c / wavelength  # angular frequency
		self.freq = c / wavelength  # frequency
		self.k = 2*np.pi / wavelength  # wavenumber
		self.energy_J = h*c / wavelength  # energy in Joules
		self.energy_ev = self.energy_J / sc.eV  # energy in electron-volts
		

class Layer:

	def __init__(self, material, num_points=0, min_wl=0, max_wl=0, thickness=0):
		self.material = material
		self.thickness = thickness
		self.num_points = num_points
		self.min_wl = min_wl  # starting wavelength
		self.max_wl = max_wl  # ending wavelength
		self.wavelength = []  # array of free space wavelengths
		self.index = []  # array of refractive indices
		self.extinct = []  # array of extinction coefficients
		self.complex = complex

	def complex_index(self, n, K):
		"""Not properly implemented"""
		self.complex = n + 1j*K


	def get_data_from_csv(self, path):
		"""Used for refractiveindex.info data"""
		with open(path, 'r') as params:
			reader = csv.reader(params)
			next(reader, None)
			for row in reader:
				wl = float(row[0])
				n = float(row[1])
				self.wavelength.append(wl)
				self.index.append(n)
				if row[2]:
					K = float(row[2])
					self.extinct.append(K)
				
	def get_data_from_txt(self, path):
		"""Used for filmetrics.com data"""
		with open(path, 'r') as params:
			header = next(params)
			unit = None
			if 'nm' in header:
				unit = 10**-3
			lines = params.readlines()
			for line in lines:
				line = line.strip().split()
				wl = float(line[0])
				if unit:
					wl = wl* unit
					print(wl)
				n = float(line[1])
				
				self.wavelength.append(wl)
				self.index.append(n)
				if line[2]:
					K = float(line[2])
					self.extinct.append(K)
				
				
	def set_wavelengths_points(self):
		new_wl = np.linspace(self.min_wl, self.max_wl, num=self.num_points, endpoint=True)
		return new_wl
		
	def interpolate_refraction(self, x0, y0):
		new_y = sp.interpolate.interp1d(x0, y0)
		return new_y
		
		
	def wavenumber(self, n, omega, theta=0):
		"""Outputs the wavenumbers for the dielectric for the given
		   angular frequency and angle"""
		k_x = n*omega/sc.c * np.cos(theta)
		k_z = n*omega/sc.c * np.sin(theta)
		return k_x, k_z

	def propagation_matrix(self, wavenumber):
		"""Inputs: wave number, thickness of medium
		   Output: propagation matrix (phase accumulation for plane wave 
					propagating through homogeneous medium)."""	
		phi = wavenumber*self.thickness
		P_i = np.matrix([[np.exp(-1j*phi), 0], [0, np.exp(1j*phi)]])
	
		return P_i

	def dynamical_matrix(self, n_, theta=0.):
		"""Inputs: index of refraction, angle of incidence
			Outputs: dynamical matrices for s-wave and p-wave."""
		
		# s-wave dynamical matrix
		m = n_ * np.cos(theta)
		Ds = np.matrix([[1, 1], [m, -m]])
					 
		# s-wave dynamical matrix
		Dp = np.matrix([[np.cos(theta), np.cos(theta)], [n_, -n_]])
		
		return Ds, Dp


def get_dict_from_yaml(yaml_file):
	"""Get data from yaml config file and put into dictionary"""

	with open(yaml_file, 'r') as yml:
		device = yaml.load(yml)
		num_layers = int(device['num_layers'])
		len_layers = int(len(device['layers']))
    
		assert_message = "num_layers ({}) and the number of layers specified ({}) in \
                          device config file do not match.".format(num_layers, len_layers)
		assert num_layers == len_layers, assert_message
	return device
    
def get_layers_from_yaml(device_dict):
	"""Takes device dictionary and outputs all layers as a list"""
	val1 = 'num_layers'
	val2 = 'num_points'
	val3 = 'min_wl'
	val4 = 'max_wl'
	
	num_layers = int(device_dict[val1])
	num_points = int(device_dict[val2])
	
	# Minimum and maximum desired wavelengths
	min_wl = float(device_dict[val3])
	max_wl = float(device_dict[val4])
	layers = []
	
	for i in range(num_layers):
	
		name = 'layer' + str(i)
		layer = device_dict['layers'][name]
		material = layer['material']
		thickness = float(layer['thickness']) * 10**-9

		layer_class = Layer(material, num_points, min_wl, max_wl, thickness)
		
		if "param_path" in layer:
			params = layer['param_path']
			if 'txt' in params:
				layer_class.get_data_from_txt(params)
			elif 'csv' in params:
				layer_class.get_data_from_csv(params)
		elif "index" in layer:
			layer_class.index = layer['index']
			layer_class.extinct = layer['extinction']
		else:
			print("Error in the yaml file.")
			
		layers.append(layer_class)

	return layers
	
	
def check_data_compatibility(layer1, layer2):
	"""NOT IMPLEMENTED.
	   Takes a list of layers obtained from get_layers_from_yaml"""
	
	row_error = "rows not all the same length in {}".format(l.material)
	column_error = "columns not all the same length in {} and {}".format(layer1.material,
																		 layer2.material)
	
	assert len(layer1.wavelength) == len(layer2.wavelength), column_error
	
	return 0

def multilayer_matrix(array):
	"""Generates the transfer matrix for a 1D array of individual
	   layer matrices"""
	TM = np.linalg.multi_dot(array)
	return TM

def inverse(matrix):
	Minv = np.linalg.inv(matrix)
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
   
def transmittance(TM, n0=0, ns=0, theta_0=0, theta_s=0):
    """Inputs: Multilayer matrix of dynamical matrices and propagation matrix."""
    M11 = TM.item((0, 0))
    t = 1/M11
    t_sq = t * np.conj(t)
    t_sq = t_sq.real

#     T = ((ns*np.cos(theta_s)) / (n0*np.cos(theta_0))) * t_sq
    T = np.linalg.det(TM) * t_sq

    return T, t


# ===== Functions below not used so much for now ==== #

def fresnel(n1, n2, k1x, k2x):
	"""Inputs:  Angular frequency of incident light
			    refractive indices for two media
			    angle of incidence (to calculate wavenumbers
	   Outputs: Reflection and transmission coefficients for s and p waves."""

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
	"""NOT IMPLEMENTED
	   Make wavelength data from desired number of points, 
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
	return wavelength, y, K
	
def filmetrics_data(data_file):
	"""Filmetrics.com T,R,A data are in nanometers"""
	wavelength = []
	Y = []
	with open(data_file) as film:
		reader = csv.reader(film, delimiter="\t")
		next(reader, None)  # Skip header
		for row in reader:
			wl = float(row[0]) * 10**-3
			wavelength.append(wl)
			Y.append(float(row[1]))
	return wavelength, Y
	

def main():

	# This stuff parses command line arguments
	parser = argparse.ArgumentParser()
	device_help = "an argument with device.yaml file in pfiles with list of layer csv files"
	parser.add_argument("device", help=device_help)
	args = parser.parse_args()
	
# 	data_Au_n = '/Users/garrek/projects/pistachio/data/Au_SiO2_n_Ciesielski.csv'
# 	data_Au_T = '/Users/garrek/projects/pistachio/data/transmittance_au_SiO2_Ciesielski.csv'
	data_Au_T = '/Users/garrek/projects/pistachio/data/Transmittance-calcs.txt.tsv'

	# Get a bunch of downloaded data
# 	wavelen_n, Au_n, Au_K = get_wave_data(data_Au_n)
# 	wavelen_T = get_wave_data(data_Au_T)[0]  # experimental data for Au transmission
# 	Au_T = get_wave_data(data_Au_T)[1]
	wavelen_T, Au_T = filmetrics_data(data_Au_T)
	
	# Inputs
	n_air = 1.0     # For testing purposes
	n_other = 1.4   # For testing purposes
	inc_ang = 0.    # If zero, p-wave and s-wave should yield same transmission
	um = 10**-6     # micrometers
	
	device = get_dict_from_yaml(args.device)  # yaml config file stored as dictionary
	layers = get_layers_from_yaml(device)  # a list of layer objects

	air = Layer('air')
	air.index = 1.0
	air.extinct = 0.
	air.complex = 1.0
	
	other = Layer('other')
	other.index = 1.4
	other.extinct = 0.
	other.complex = 1.4
	
	au = layers[0]
	Au_n = au.index
	Au_K = au.extinct
	new_wl = au.set_wavelengths_points()
	au.index = au.interpolate_refraction(au.wavelength, au.index)(new_wl)
	au.extinct = au.interpolate_refraction(au.wavelength, au.extinct)(new_wl)
	au.wavelength = new_wl
		
	print("smallest wavelength", au.wavelength[0])
	print("largest wavelength", au.wavelength[-1])

	# Outputs
	wavelen = [i for i in layers[0].wavelength]
	print("num of wavelengths", len(wavelen))

	print('index points', len(layers[0].index))
	R = []  # calculated reflectance data
	T = []  # calculated transmittance data

	num_points = device['num_points']   # TODO: This is really sketchy.
	num_layers = device['num_layers']
	print('number of layers:', num_layers)
	print("number of points", num_points)

	for i in range(num_points):
		#TODO: Deal with data that have a different number of points 
		M_list = []  # list of matrices to be multiplied
		lmbda = wavelen[i] * um
		light = Light(lmbda)  # Make an instance of Light class
		omega = light.omega

		# Add air to matrix list (reflection region)		
		air.wavelength = lmbda
		D0 = air.dynamical_matrix(air.index)[0]
		D0inv = inverse(D0)
		M_list.append(D0inv)

		for layer in layers:
			#Cycle through layer object instances stored in the 'layers' list above.
			
			n = layer.index[i] + 1j*layer.extinct[i]
			kx = layer.wavenumber(n, omega)[0]
			
			D = layer.dynamical_matrix(n)[0]
			Dinv = inverse(D)
			P = layer.propagation_matrix(kx)
					
			M_list.extend([D, P, Dinv])

		# Add other to matrix list (transmission region)
		other.wavelength = lmbda
		Ds = other.dynamical_matrix(other.index)[0]
		M_list.append(Ds)

# 		M_i = multilayer_matrix(M_list)
		M_i = multilayer_matrix([D0inv, D, P, Dinv, Ds])
					
		trn = transmittance(M_i, air.index, other.index)[0].real
		ref = reflectance(M_i)[0]
		T.append(trn)
		R.append(ref)

		
	
	# Make Plots
	fig, axs = plt.subplots(3, 1, sharex=True)
	
	ax = axs[0]
# 	ax.plot(wavelen, Au_n, label="n, downloaded data")
	ax.plot(au.wavelength, au.index, label="n, downloaded")
	ax.plot(au.wavelength, au.extinct, label="K, downloaded data")
	ax.set_ylabel('refractive index for Au')
	ax.legend()
	
	ax = axs[1]
	ax.plot(au.wavelength, T, label="calculated data")
	ax.plot(wavelen_T, Au_T, label="downloaded data")
	ax.set_ylabel('transmittance')
# 	ax.set_xlim(0, 5)
	ax.legend()
	
	ax = axs[2]
	ax.plot(au.wavelength, R, label="calculated data")
	ax.set_xlabel('wavelength ($\mu$m)')
	ax.set_ylabel('reflectance')
	ax.legend()
	
# 	ax = axs[3]
# 	ax.plot(au_x, au_y(au_x), label="testing interpolation")
# 	ax.legend()
	
	fig.suptitle('Index of refraction, transmission, reflection for {}'.format(au.material))
	fig.tight_layout()
	plt.show()
	

if __name__ == '__main__':
	main()
