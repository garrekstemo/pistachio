#! /anaconda3/bin/python

# Convention used
# Psi(x, t) = Psi_0 * exp(i(kx - wt))
# n = n' + i*n'', where n' is real part of refractive index and n'' is imaginary part.

import argparse
import sys
import csv
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.constants as sc
import scipy.interpolate
from TMM_tests import *
from ruamel_yaml import YAML
# import pdb

c = sc.c  # speed of light
h = sc.h  # planck's constant
yaml = YAML()


class Light:

	def __init__(self, wavelength):
		
		self.wavelength = wavelength
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
		
	def __repr__(self):
		a = "{} \n".format(self.material)
		b = "thickness: {}\n".format(self.thickness)
		c = "wavelength: {}\n".format(self.wavelength)
		d = "refractive index: {}\n".format(self.index)
		e = "extinction coefficient: {}\n".format(self.extinct)
		return a+b+c+d+e
				

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
				n = float(line[1])
				
				self.wavelength.append(wl)
				self.index.append(n)
				if line[2]:
					K = float(line[2])
					self.extinct.append(K)
				
	def set_wavelengths(self):
		new_wl = np.linspace(self.min_wl, self.max_wl, num=self.num_points, endpoint=True)
		return new_wl
		
	def interpolate_refraction(self, x0, y0):
		# WARNING: fill_value='extrapolate' might cause problems!!
		new_y = sp.interpolate.interp1d(x0, y0, fill_value="extrapolate")
		return new_y

	def make_new_data_points(self):
		"""Makes new data points based on user-defined num_points and interpolation."""
		if not isinstance(self.index, list):
			self.index = [self.index]*self.num_points	
			self.extinct = [self.extinct]*self.num_points
			self.wavelength = self.set_wavelengths()
		elif isinstance(self.index, list):
			new_n = self.interpolate_refraction(self.wavelength, self.index)
			new_K = self.interpolate_refraction(self.wavelength, self.extinct)
			self.wavelength = self.set_wavelengths()
			self.index = new_n(self.wavelength)
			self.extinct = new_K(self.wavelength)

		
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
		
		if args.swave:
			return Ds
		elif args.pwave:
			return Dp
		else:
			return "No wave polarization passed in."
		


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
			layer_class.wavelength = layer['wavelength']
		else:
			print("Error in the yaml file.")
			
		layers.append(layer_class)

	return layers


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

def wavelengths_loop(pts, wav_list, substrate, layers, air, inc_ang):
	
	T = []
	R = []
	M = []
	
	um = 10**-6
	
	for i in range(pts):
		
		M_list = []
		lmbda = wav_list[i] * um
		light = Light(lmbda)
		omega = light.omega
		
		# Add substrate to matrix list (reflection region)
		n = substrate.index[i] + 1j*substrate.extinct[i]
		D0 = substrate.dynamical_matrix(n, inc_ang)
		D0inv = inverse(D0)
		M_list.append(D0inv)
		
		for layer in layers:
			n  = layer.index[i] + 1j*layer.extinct[i]
			kx = layer.wavenumber(n, omega, inc_ang)[0]
			
			D = layer.dynamical_matrix(n, inc_ang)
			Dinv = inverse(D)
			P = layer.propagation_matrix(kx)
			
			M_list.extend([D, P, Dinv])
		
		# Add transmission region medium to matrix list
		Ds = air.dynamical_matrix(air.index, inc_ang)
		M_list.append(Ds)
		M_i = multilayer_matrix(M_list)
		
		trn = transmittance(M_i, substrate.index[i], air.index)[0]
		ref = reflectance(M_i)[0]
		T.append(trn)
		R.append(ref)
		M.append(M_i)
	
	return M, T, R

def new_em_coeff(M, A0, B0):
	"""Takes transfer matrix, initial magnitude of EM wave, A0, B0.
	Outputs transformed magnitudes, As' and Bs'"""
	
	E0 = np.array([A0, B0])
	Es = inverse(M).dot(E0)
	return Es
	

# ===== Fresnel equation functions below not used so much for now ==== #

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
	
def reference_data(data_file):
	"""Gets reference data downloaded from
	websites. Filmetrics.com T,R,A data are in nanometers"""
	wavelength = []
	Y = []
	unit = 1  # sets order of magnitude (nm or um)
	with open(data_file, 'r') as ref:
	
		reader = None
		if 'txt.tsv' in str(ref):
			print('Filmetrics data. Units in nm')
			unit = 10**-3
			reader = csv.reader(ref, delimiter="\t")
		else:
			print("Refractiveindex.info data. Units in um")
			reader = csv.reader(ref)
		next(reader, None)  # Skip header
		for row in reader:
			wl = float(row[0]) * unit  # MAKE SURE UNITS ARE CORRECT
			wavelength.append(wl)
			Y.append(float(row[1]))
	return wavelength, Y


def main():

	print("\nStart simulation\n")

	# Store a bunch of downloaded data
	if args.trns:
		wl_T_data, T_data = reference_data(args.trns)
	if args.refl:
		wl_R_data, R_data = reference_data(args.refl)
	else:
		print("Using no reference data.\n")
		
	# Inputs
	um = 10**-6     # micrometers
	
	device = get_dict_from_yaml(args.device)  # yaml config file stored as dictionary
	layers = get_layers_from_yaml(device)  # a list of layer objects
	
	# If zero, p-wave and s-wave should yield same transmission
	em_wave = device['wave']
	inc_ang = em_wave['theta_in'] #np.pi/3

	# TODO: put this in yaml config file
	air = Layer('air')
	air.index = 1.0
	air.extinct = 0.
	air.complex = 1.0
	
	# Get bounding values and populate each layer with data
	num_points = device['num_points']
	num_layers = device['num_layers']
	MIN_WL = device['min_wl']
	MAX_WL = device['max_wl']
	# Now we know how many points there are and we can generate some wavelengths
	# based on min and max wavelengths in our config file.

	print('number of layers:', num_layers)
	print("number of points", num_points)
	print("Starting wavelength (um)", MIN_WL)
	print("Ending wavelength (um)", MAX_WL)
	
	sub_name = device['substrate']['material']
	substrate = Layer(sub_name, num_points, MIN_WL, MAX_WL)

	if 'param_path' in device['substrate']:
		sub_params = device['substrate']['param_path']
		substrate.get_data_from_csv(device['substrate']['param_path'])
		substrate.make_new_data_points()
	else:
		substrate.index = device['substrate']['index']
		substrate.extinct = device['substrate']['extinction']
		substrate.make_new_data_points()
	
	for layer in layers:
		# interpolating from downloaded index data, making new data points
		print("\n", str(layer.material) + ", d=" + str(int(layer.thickness*10**9)) + "nm")
		layer.make_new_data_points()

	wavelen = layers[0].wavelength  # still in units of um
	
	# Outputs
	
	M, T, R = wavelengths_loop(num_points, wavelen, substrate, layers, air, inc_ang)
	
	for m in M:
		new_em_coeff(m, em_wave['A0'], em_wave['B0'])
	
	# Make Plots
	print("_"*50)
	print("Generating plots...")
	fig, axs = plt.subplots(2, 1, sharex=True)
	
	ax = axs[0]
	if args.trns:
		ax.plot(wl_T_data, T_data, color='0.3', label="downloaded data")
	ax.plot(wavelen, T, 'r--', label="calculated data")
	ax.set_ylabel('transmittance')
# 	ax.set_xlim(4, 15.)
	ax.legend()
	
	ax = axs[1]
	ax.plot(wavelen, R, 'r-', label="calculated data")
	if args.refl:
		ax.plot(wl_R_data, R_data, '--', label="downloaded data")
	ax.set_xlabel('wavelength ($\mu$m)')
	ax.set_ylabel('reflectance')
	ax.legend()
	
# 	d = str(int(au.thickness*10**9)) + ' nm'
	title = "Index of refraction (n, K), transmission, reflection for 1000nm SiO2 and 35nm Au" #for {} {}.format(d, au.material)
# 	FPI_title = "Transmission, reflection for SiO2 (100nm), Au (35nm), Air (5um), Au, SiO2"
	plt.suptitle(title)
	plt.tight_layout()
	plt.subplots_adjust(top=0.9)
	plt.show()
	

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	device_help = "path for a yaml file from config_files describing a device"
	reflectance_help = "path for reflectance data downloaded from filmetrics"
	transmittance_help = "path for transmittance data downloaded from filmetrics"
	pwave_help = "boolean, calculate for p-wave"
	swave_help = "boolean, calculate for s-wave"
	
	parser.add_argument('--debug',
						action='store_true',
						help="doesn't do anything right now")
	parser.add_argument("device",
						help=device_help)
	parser.add_argument('-p', '--pwave', help=pwave_help, action='store_true')
	parser.add_argument('-s', '--swave', help=swave_help, action='store_true')
	parser.add_argument('-T', '--trns',
						help=transmittance_help)
	parser.add_argument('-R', '--refl',
						help=reflectance_help)
	debug_mode = parser.parse_args().debug
	if debug_mode:
		print("Debug mode\n")
	args = parser.parse_args()
	
	main()
