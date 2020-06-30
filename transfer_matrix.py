#! /anaconda3/bin/python

# Convention used
# Psi(x, t) = Psi_0 * exp(i(kx - wt))
# n = n' + i*n'', where n' is real part of refractive index and n'' is imaginary part.

import argparse
import os
import sys
import csv
from itertools import tee
import logging
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy as sp
import scipy.constants as sc
import scipy.interpolate
import time
from ruamel_yaml import YAML
import pdb


c = sc.c  # speed of light
h = sc.h  # planck's constant
yaml = YAML()


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
		self.waves = {}  # Needs to be curly braces
		
	def __repr__(self):
		a = "{} \n".format(self.material)
		b = "thickness: {}\n".format(self.thickness)
		c = "wavelength: {}\n".format(self.wavelength)
		d = "refractive index: {}\n".format(self.index)
		e = "extinction coefficient: {}\n".format(self.extinct)
		return a+b+c+d+e

	def get_data_from_csv(self, index_path):
		"""Used for refractiveindex.info data. This site uses um for wavelength units."""
		with open(index_path, 'r') as params:
			reader = csv.reader(params)
			next(reader, None)
			for row in reader:
				wl = float(row[0])
				n = float(row[1])
				self.wavelength.append(wl)
				self.index.append(n)
				try:
					K = float(row[2])
					self.extinct.append(K)
				except IndexError:
					self.extinct.append(0.)
				
	def get_data_from_txt(self, index_path):
		"""Used for filmetrics.com data"""
		with open(index_path, 'r') as params:
			header = next(params)
			unit = 1
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
		"""Creates a new set of linearly-spaced wavelengths from start and 
		   endpoint wavelengths and number of points specified in yaml config file."""
		new_wl = np.linspace(self.min_wl, self.max_wl, num=self.num_points, endpoint=True)
		return new_wl
		
	def interpolate_refraction(self, x0, y0):
		#WARNING: fill_value='extrapolate' might cause problems!!
		new_y = sp.interpolate.interp1d(x0, y0, fill_value="extrapolate")
		return new_y

	def make_new_data_points(self):
		"""Makes new data points based on user-defined num_points and interpolation."""
		if not isinstance(self.index, list):
			self.index = [self.index]*self.num_points	
			self.extinct = [self.extinct]*self.num_points
			self.wavelength = self.set_wavelengths()
			for idx, lmbda in enumerate(self.wavelength):
				self.waves[lmbda] = self.index[idx] + self.extinct[idx]

		elif isinstance(self.index, list):
			new_n = self.interpolate_refraction(self.wavelength, self.index)
			new_K = self.interpolate_refraction(self.wavelength, self.extinct)
			self.wavelength = self.set_wavelengths()
			self.index = new_n(self.wavelength)
			self.extinct = new_K(self.wavelength)
			for idx, lmbda in enumerate(self.wavelength):
				self.waves[lmbda] = self.index[idx] + 1j*self.extinct[idx]
		
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
		P_i = np.array([[np.exp(-1j*phi), 0], [0, np.exp(1j*phi)]])
	
		return P_i

	def dynamical_matrix(self, n_, theta=0.):
		"""Inputs: index of refraction, angle of incidence
			Outputs: dynamical matrices for s-wave and p-wave."""
		
		# s-wave dynamical matrix
		m = n_ * np.cos(theta)
		Ds = np.array([[1, 1], [m, -m]])
					 
		# p-wave dynamical matrix
		Dp = np.array([[np.cos(theta), np.cos(theta)], [n_, -n_]])
		
		if args.swave:
			return Ds
		elif args.pwave:
			return Dp
		else:
			print("")
			print("Error: No wave polarization passed in. See --help for more info. Exiting....")
			sys.exit()


def get_dict_from_yaml(yaml_file):
	"""Get data from yaml config file and put into dictionary"""

	with open(yaml_file, 'r') as yml:
		device = yaml.load(yml)
	return device

def get_layers_from_yaml(device_dict):
	"""Takes device dictionary and outputs all layer objects as a list."""
	val1 = 'layers'
	val2 = 'num_points'
	val3 = 'min_wl'
	val4 = 'max_wl'
	
	num_layers = len(device_dict[val1])
	print("Number of layers:", num_layers)
	num_points = int(device_dict[val2])
	
	# Minimum and maximum desired wavelengths
	min_wl = float(device_dict[val3])
	max_wl = float(device_dict[val4])
	layers = []
	
	for i in range(num_layers):
	
		layer_str = 'layer' + str(i)
		layer = device_dict['layers'][layer_str]
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
			print("ERROR: Incorrect yaml config format. Reference default template.")
			
		layers.append(layer_class)

	return layers


def get_beam_profile(beam_csv):
	"""Gets field distribution data from FTIR csv file.
	   Outputs list of wavenumbers and field amplitudes.
	   WARNING: Jasco FTIR has weird non utf-8 characters in the footer and 
	   partway through data. This annoyance is kind of dealt with. Might have problems."""
	
	with open(beam_csv, 'r', encoding="utf8", errors="ignore") as beam:
		num_header_rows = 18
		num = 0
		reader = csv.reader(beam)
		wavenum = []
		field = []
		
		while num <= num_header_rows:
			next(reader, None)
			num += 1

		next_reader, reader = tee(reader)
		next(next_reader)
		try:
			for row in reader:
				try:
					x = float(row[0])
					y = float(row[1])
					wavenum.append(x)
					field.append(float(y))
					next(next_reader)
				except (IndexError, ValueError):
					# This is where we deal with weird utf8 characters
					pass
		except StopIteration:
			pass
		
	return wavenum, field

def attenuation_coefficient(beam_intensity, x):
	"""NOT PROPERLY IMPLEMENTED. 
	Takes beam intensity as array"""
	
	dI = np.gradient(beam_intensity, x)
	
	return dI / beam_intensity
	
def kramers_kronig(alpha):
	"""NOT IMPLEMENTED.
	Takes attenuation coefficient. Returns real part of index of refraction."""

def reflectance(M_):
    """Input: multilayer matrix, M.
       Output: reflectance calculation."""
    M21 = M_.item((1, 0))
    M11 = M_.item((0, 0))
    r = M21 / M11
    r_sq = r * np.conj(r)
    r_sq = r_sq.real
    R = r_sq

    return R
   
def transmittance(TM):
    """Inputs: Multilayer matrix of dynamical matrices and propagation matrix."""
    M11 = TM.item((0, 0))
    t = 1/M11
    t_sq = t * np.conj(t)
    t_sq = t_sq.real

    T = np.linalg.det(TM) * t_sq

    return T

def build_matrix_list(wavelength, theta, layers):
	"""
	Makes a wavelength object and bounds. Outputs a list of matrices that will
	make the transfer matrix for that wavelength of light propagating 
	through the device."""
	
	compute_wl = wavelength * 10**(-6)
	omega = 2 * np.pi * sc.c / compute_wl
	matrices = []

	for idx, layer in enumerate(layers):

		n = layer.waves[wavelength]
		kx = layer.wavenumber(n, omega, theta)[0]
		D = layer.dynamical_matrix(n, theta)
		Dinv = np.linalg.inv(D)
		P = layer.propagation_matrix(kx)
		
		if idx == 0:
			matrices.append(Dinv)
			
		elif idx == len(layers)-1:
			matrices.append(D)

		else:
			matrices.append(D)
			matrices.append(P)
			matrices.append(Dinv)
			
	return matrices


def matrix_product(matrices):
	"""Product of matrices"""
	
	i = 0
	M = np.array([[1,0], [0,1]])
	for i in matrices:
		M = np.dot(M, i)
	return M

def field_amp(matrix_list, A0_, B0_):
	"""E_s = M*E_0
		= D0inv * (Di_inv*P*Di)^n * D0  * E_0
	Takes transfer matrix, final magnitude of EM wave, As, Bs.
	Outputs E-field amplitude for given matrix. This approach starts calculations from
	the incident wave."""
	
	field_rev = []  #reversed order field amplitudes
	E0 = np.array([A0_, B0_])

	i = len(matrix_list) - 1
	slice = 3  # slice first three matrices from list after each dot product
	M = matrix_product(matrix_list)
	E_s = np.dot(M, E0)
	field_rev.append(E_s)
	matrix_list.pop(0)
	
	while i != 1:
		M_i = matrix_product(matrix_list)
		E_i = np.dot(M_i, E0)
		field_rev.append(E_i)
		matrix_list = matrix_list[slice:]
		i -= slice
	if len(matrix_list) == 1:
		E_i = np.dot(matrix_list[0], E0)
		field_rev.append(E_i)

	field = list(reversed(field_rev))

	return field
	

def output_TRA(wavelens, trans, refl, absor):
	"""Writes transmission, reflectance, abosorbance data to csv file"""
	
	print("")
	output = os.path.abspath(args.output)
	with open (output, 'w') as out_file:
	#TODO: Also output wavenumbers
		filewriter = csv.writer(out_file, delimiter=',')
		header = ['Wavelength', 
				  'Transmittance', 
				  'Reflectance', 
				  'Absorptance']
		filewriter.writerow(header)
		i = 0
		while i < len(trans):
			row = [wavelens[i], trans[i], refl[i], absor[i]]
			filewriter.writerow(row)
			i+=1
	print("Wrote results to {}".format(output))
	return 0
	
	
def output_field_profile(wavelens, layers, E_amps):
	"""Take list of wavelengths"""
	k = 0.2  # test wavenumber, 2000 cm^-1
	a = E_amps[100]
	layer_coords = []
	for idx, layer in enumerate(layers):
		if idx == 0:
			layer_coords.append(0)
		else:
			d = int(layer.thickness * 10**9)
			x_i = layer_coords[idx-1] + d
			layer_coords.append(x_i)

	# Generate x-coordinates
	num_pts = 10000
	MAX = layer_coords[-1] + 0.1*layer_coords[-1]
	LAST = len(layer_coords) - 1
	x = []
	profile = []
	for idx, l in enumerate(layer_coords):
		if idx == LAST:
# 			print(idx, "last one")
			domain = np.linspace(l, MAX, num_pts)
			x.append(domain)
		else:
			domain = np.linspace(l, layer_coords[idx+1], num_pts)
			x.append(domain)

	fig, ax = plt.subplots()
	
	# Generate E-field profile
	for idx, x_i in enumerate(x):
		f = a[idx] * np.exp(1j * k * x_i)
		profile.append(f)
		ax.plot(x_i, f, label='{}'.format(layers[idx].material))

	ax.legend()

	ax.axvline(x=layer_coords[1])
	ax.axvline(x=layer_coords[2])
	ax.axvline(x=layer_coords[3])


	field_output = '/Users/garrek/projects/pistachio/data/out/field_test.csv'
	with open (field_output, 'w') as out_file:
		filewriter = csv.writer(out_file, delimiter=',')
		header = ['x', 'field']
		filewriter.writerow(header)
		i = 0
		while i < len(x):
			row = [x[i], profile[i]]
			filewriter.writerow(row)
			i+=1
	print("Wrote field output to {}".format(field_output))		

	plt.show()
	
	return 0

# ===== Fresnel equation functions below not used for now ===== #
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
	D_ij = 1/t_ij * np.array([[1, r_ij], [r_ij, 1]])
	return D_ij
# ========= ========= ========= ========= ========== ========= ======== #


def main_loop():
	"""Executes transfer matrix and other functions
	   Writes output file."""

	# Inputs
	device = get_dict_from_yaml(args.device)  # yaml config file stored as dictionary
	layers = get_layers_from_yaml(device)  # a list of layer objects
	em_wave = device['wave']
	inc_ang = em_wave['theta_in']
	
	for layer in layers:
		# interpolating from downloaded index data, making new data points
		#TODO: Why am I doing this?
		print(str(layer.material) + ", d=" + str(int(layer.thickness*10**9)) + "nm")
		layer.make_new_data_points()

	# We should separate this from the Layer class.
	# This should be stored as an array somewhere.
	wavelens = layers[0].wavelength  # still in units of um

	
	#TODO: Organize outputs better cuz there'll probably be a lot of them.
	# Outputs
	T = []
	R = []
	A = []
	E_amps = []
	intensity = []
	
# 	alpha = attenuation_coefficient(efield, k)

	for lmbda in wavelens:

		M = build_matrix_list(lmbda, inc_ang, layers)
		TM = np.linalg.multi_dot(M)
		trns = transmittance(TM).real
		refl = reflectance(TM).real
# 		field = field_amp(M, device["wave"]['A0'], device['wave']['B0'])
		logging.info('Only using forward-propagating field values')
# 		field = [f[0] for f in field]
			
# 		E_amps.append(field)
		T.append(trns)
		R.append(refl)
		abso = 1 - trns - refl
		A.append(abso)

	#Write everything to a csv file
	output_TRA(wavelens, T, R, A)
# 	output_field_profile(wavelens, layers, E_amps)


def main():

	LOG_FILENAME = 'transfer_matrix.log'
	logging.basicConfig(filename=LOG_FILENAME, level=logging.INFO)
	print("\nStart simulation\n")
	
	start_time = time.time()
	main_loop()
	end_time = time.time()
	elapsed_time = np.round(end_time - start_time, 4)
	logging.info('elapsed time: {} seconds'.format(elapsed_time))
	
	print("Simulation time: {} seconds".format(elapsed_time))


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	device_help = "path for a yaml file from config_files describing a device"
	output_help = "path and name for output data file (must be .csv)"
	pwave_help = "boolean, calculate for p-wave"
	swave_help = "boolean, calculate for s-wave"
	
	parser.add_argument('--debug',
						action='store_true',
						help="doesn't do anything right now")
	parser.add_argument("device",
						help=device_help)
	parser.add_argument("output", help=output_help)
	parser.add_argument('-p', '--pwave', help=pwave_help, action='store_true')
	parser.add_argument('-s', '--swave', help=swave_help, action='store_true')
# 	parser.add_argument('-T', '--angles', action='store_true')
	debug_mode = parser.parse_args().debug
	if debug_mode:
		print("Debug mode\n")
		
	args = parser.parse_args()
	
	output_message = "Output file must be .csv"
	assert args.output[-4:] == ".csv", output_message
		
	main()
