#!/usr/bin/env python

"""
Name: Transfer Matrix
Author: Garrek Stemo
Affiliation: Nara Institute of Science and Technology


Convention used:
Psi(x, t) = Psi_0 * exp(i(kx - wt))
n = n' + i*n'', where n' is real part of refractive index and n'' is imaginary part.
Yeh, Pochi. 2005. Optical Waves in Layered Media.
"""

import argparse
import codecs
import csv
import importlib.resources as pkg_resources
import logging
import os
import pdb
import sys
import time
from itertools import tee
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import multiprocessing
import numpy as np
from ruamel_yaml import YAML
import scipy as sp
import scipy.constants as sc
import scipy.interpolate
from tqdm import tqdm
import data.refractive_index_data  # import directory containing refractive index info
import results

yaml = YAML()
# FORMATTER = logging.Formatter("%(asctime)s — %(name)s — %(levelname)s - %(message)s")
FORMATTER = logging.Formatter("%(message)s")
LOG_FILE = 'transfer_matrix.log'

class Wave:
	"""
	Contains information about incident light and includes functions to
	perform calculations and unit conversion.
	"""
	def __init__(self, min_wl, max_wl, num_points):
		self.min_wl = min_wl  # Starting wavelength
		self.max_wl = max_wl  # Ending wavelength
		self.num_points = num_points  # Number of wavelengths between min_wl and max_wl
		self.wavelengths = []

	def make_wavelengths(self):
		"""Generate linearly-spaced wavelengths between min_wl and max_wl
		   with number of wavelengths equal to num_points"""
		self.wavelengths = np.linspace(self.min_wl, self.max_wl, self.num_points)

	def wavelength_to_wavenumber(self, wavelengths):
		"""Convert micrometers to cm-1 for an array of wavelengths"""
		return [10**4 / l for l in wavelengths]

	def wavenumber_to_eV(self, wavenumbers):
		"""Convert from cm-1 to eV for an array of wavenumbers"""
		return [sc.h * sc.c * k * 100.0 / sc.elementary_charge for k in wavenumbers]


class Layer:
	"""
	Contains information and functions related to a single material layer in
	a multi-layer device.
	"""
	def __init__(self, material, num_points, min_wl=0.0, max_wl=0.0, thickness=0.0):
		self.material = material			# Material name string
		self.thickness = thickness			# Layer thickness (be consistent, but test files are in nm)
		self.num_points = num_points		# Number of wavelengths to scan through
		self.wavelengths = []  				# wavelengths from refractive index data. Used only for testing.
		self.refractive_index = [] 			# Array of refractive indices (real part)
		self.extinction_coeff = []  		# Array of extinction coefficients (imaginary part)
		self.complex_refractive = {}  		# Needs to be curly braces

	def __repr__(self):
		a = "{} \n".format(self.material)
		b = "thickness: {}\n".format(self.thickness)
		c = "wavelengths: {}\n".format(self.wavelengths)
		d = "refractive index: {}\n".format(self.refractive_index)
		e = "extinction coefficient: {}\n".format(self.extinction_coeff)
		return a+b+c+d+e

	def get_data_from_csv(self, refractive_filename):
		"""
		Extract refractive index data from file downloaded from refractiveindex.info. 
		This site uses micrometers for wavelength units.
		"""

		with pkg_resources.path(data.refractive_index_data, refractive_filename) as params:
			# pkg_resources will return a path, which includes the csv file we want to read
			params = os.path.abspath(params)
			with open(params, 'r', encoding='utf-8') as csv_file:
				csvreader = csv.reader(csv_file)
				next(csvreader, None)
				for row in csvreader:
					wl = float(row[0])
					n = float(row[1])
					self.wavelengths.append(wl)
					self.refractive_index.append(n)
					try:
						K = float(row[2])
						self.extinction_coeff.append(K)
					except IndexError:
						self.extinction_coeff.append(0.0)

	def get_data_from_txt(self, index_path):
		"""
		Extract refractive index data from file downloaded from for filmetrics.com.
		WARNING: Not tested in a long time. The format may have changed. Code may have changed.
		"""
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

				self.wavelengths.append(wl)
				self.refractive_index.append(n)
				if line[2]:
					K = float(line[2])
					self.extinction_coeff.append(K)

	def make_new_data_points(self, wavelengths):
		"""
		Makes new data points based on user-defined num_points and interpolation.
		Input: list of wavelengths in micrometers generated from yaml config file.
		Output: refractive index values (real and imaginary) mapped to wavelengths.
		"""
		
		if not isinstance(self.refractive_index, list):
			self.refractive_index = [self.refractive_index]*self.num_points
			self.extinction_coeff = [self.extinction_coeff]*self.num_points
			for idx, lmbda in enumerate(wavelengths):
				self.complex_refractive[lmbda] = self.refractive_index[idx] + self.extinction_coeff[idx]

		elif isinstance(self.refractive_index, list):
			new_n = sp.interpolate.interp1d(self.wavelengths, self.refractive_index, fill_value='extrapolate')
			new_K = sp.interpolate.interp1d(self.wavelengths, self.extinction_coeff, fill_value='extrapolate')
			self.refractive_index = new_n(wavelengths)
			self.extinction_coeff = new_K(wavelengths)
			for idx, lmbda in enumerate(wavelengths):
				self.complex_refractive[lmbda] = self.refractive_index[idx] + 1j*self.extinction_coeff[idx]

	def get_wavenumber(self, n, omega, theta=0):
		"""Outputs the wavenumbers for the dielectric for the given
		   angular frequency and angle"""

		k_x = n*omega/sc.c * np.cos(theta)
		k_z = n*omega/sc.c * np.sin(theta)
		return k_x, k_z


def get_dict_from_yaml(yaml_file):
	"""Get data from yaml config file and put into dictionary"""

	with open(yaml_file, 'r') as yml:
		device = yaml.load(yml)
	return device

def get_layers_from_yaml(device_dict):
	"""
	Takes device dictionary from yaml and outputs all layer objects as a list.
	Input: dictionary with multi-layer device data already obtained rom yaml.load()
	Output: List of Layer classes with parameters pulled from yaml config file.
	"""
	key1 = 'layers'
	key2 = 'num_points'
	key3 = 'min_wavelength'
	key4 = 'max_wavelength'
	key5 = 'wave'

	num_layers = len(device_dict[key1])
	num_points = int(device_dict[key2])

	# Minimum and maximum wavelengths from yaml config file
	min_wl = float(device_dict[key3])
	max_wl = float(device_dict[key4])
	theta_i = device_dict[key5]['theta_i']
	theta_f = device_dict[key5]['theta_f']
	num_angles = device_dict[key5]['num_angles']

	print("Wavelengths: {} in range [{}, {}]".format(num_points, min_wl, max_wl))
	print("Sweep through {} angles in range [{}, {}]".format(num_angles, theta_i, theta_f))
	print('_'*50)
	print('Device Configuration')
	layers = []

	for i in range(num_layers):

		layer_str = 'layer' + str(i)
		layer = device_dict['layers'][layer_str]
		material = layer['material']
		thickness = float(layer['thickness']) * 10**-9
		layer_class = Layer(material, num_points, min_wl, max_wl, thickness)

		if "refractive_filename" in layer:
			params = layer['refractive_filename']
			if 'txt' in params:
				layer_class.get_data_from_txt(params)
			elif 'csv' in params:
				layer_class.get_data_from_csv(params)
		elif "refractive_index" in layer:
			layer_class.refractive_index = layer['refractive_index']
			layer_class.extinction_coeff = layer['extinction_coeff']
			layer_class.wavelengths = layer['wavelength']
		else:
			print("ERROR: Incorrect yaml config format. Reference default template.")

		layers.append(layer_class)
		print(str(layer_class.material) + ", d=" + str(int(layer_class.thickness*10**9)) + "nm")
	print('_'*50)
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


def propagation_matrix(wavenumber, layer_thickness):
	"""Inputs: wave number, thickness of medium
	   Output: propagation matrix (phase accumulation for plane wave
				propagating through homogeneous medium)."""
	phi = wavenumber*layer_thickness
	P_i = np.array([[np.exp(-1j*phi), 0.0], [0.0, np.exp(1j*phi)]])

	return P_i


def dynamical_matrix(n_, theta=0.0, wave_type='mixed'):
	"""Inputs: index of refraction, angle of incidence
		Outputs: dynamical matrices for s-wave and p-wave."""

	# s-wave dynamical matrix
	m = n_ * np.cos(theta)
	Ds = np.array([[1, 1], [m, -m]])

	# p-wave dynamical matrix
	Dp = np.array([[np.cos(theta), np.cos(theta)], [n_, -n_]])

	if wave_type == 's-wave':
		return Ds
	elif wave_type == 'p-wave':
		return Dp
	elif wave_type == 'mixed':
		pol_msg = "\nNo acceptable polarization passed in.\nMixed-wave not yet supported.\nSee --help for more info. Exiting...."
		print(pol_msg)
		sys.exit()


def find_reflectance(M_):
    """Input: multilayer matrix, M.
       Output: Reflectance calculated from elements of transfer matrix."""
    M21 = M_.item((1, 0))
    M11 = M_.item((0, 0))
    r = M21 / M11
    r_sq = r * np.conj(r)
    r_sq = r_sq.real
    R = r_sq

    return R


def find_transmittance(M_):
    """
    Inputs: Multilayer matrix of dynamical matrices and propagation matrix.
    Output: Transmittance calculated from elements of transfer matrix.
    """
    M11 = M_.item((0, 0))
    t = 1/M11
    t_sq = t * np.conj(t)
    t_sq = t_sq.real

    T = np.linalg.det(M_) * t_sq

    return T


def build_matrix_list(wavelength, theta, layers, wave_type):
	"""
	Makes a wavelength object and bounds. Outputs a list of matrices that will
	make the transfer matrix for that wavelength of light propagating
	through the device."""

	compute_wl = wavelength * 10**(-6)
	omega = 2 * np.pi * sc.c / compute_wl
	matrices = []

	for idx, layer in enumerate(layers):

		n = layer.complex_refractive[wavelength]
		kx = layer.get_wavenumber(n, omega, theta)[0]
		D = dynamical_matrix(n, theta, wave_type)
		Dinv = np.linalg.inv(D)
		P = propagation_matrix(kx, layer.thickness)

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
		M = np.matmul(M, i)
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
	slice_matrices = 3  # slice first three matrices from list after each dot product

	M = matrix_product(matrix_list)
	E_s = np.dot(M, E0)
	field_rev.append(E_s)
	matrix_list.pop(0)

	while i != 1:
		M_i = matrix_product(matrix_list)
		E_i = np.dot(M_i, E0)
		field_rev.append(E_i)
		matrix_list = matrix_list[slice_matrices:]
		i -= slice_matrices

	if len(matrix_list) == 1:
		E_i = np.dot(matrix_list[0], E0)
		field_rev.append(E_i)

	field = list(reversed(field_rev))

	return field


def write_tmm_results(angle, output_dir, rows):
	"""Writes transmission, reflectance, abosorbance data to csv file"""

	wavelens, trans, refl, absor = rows
	file_name = 'deg' + str(angle) + '_.csv'  # Underscore is to make it work with reading in angle in another method
	output_file = os.path.join(output_dir, file_name)

	with open (output_file, 'w', encoding='utf8', newline='') as out_file:

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


def output_field_profile(wavelens, layers, E_amps):
	"""
	Take list of wavelengths
	WARNING: Needs to be tested.
	"""
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

	field_output = ''
	with pkg_resources.path('data', 'out/field_test.csv') as field:
		field_output = field
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


def perform_transfer_matrix(sim_path, angle, wavelengths, layers, wave_type):

	transmittance = []
	reflectance = []
	absorbance = []

	for lmbda in wavelengths:
		theta = angle * np.pi / 180.0
		M = build_matrix_list(lmbda, theta, layers, wave_type)
		TM = np.linalg.multi_dot(M)
		T = find_transmittance(TM).real
		R = find_reflectance(TM).real

		transmittance.append(T)
		reflectance.append(R)
		absorbance.append(1 - T - R)

	#Write everything to a csv file
	row = [wavelengths, transmittance, reflectance, absorbance]
	write_tmm_results(angle, sim_path, row)
# 	results = [wavelengths, transmittance, reflectance, absorbance]
# 	return angle, sim_path, results


def angle_resolved_multiprocess(device_yaml, output_dir, wave_type):
	"""
	Inputs: yaml file containing information about device and incident radiation.
	Outputs: Executes transfer matrix and other functions
	   		 Writes output file.
	   		 
	This process uses the Python multiprocessing library. Each angle is assigned
	its own process to speed up transfer matrix calculations for angle-tuned
	simulations.
	"""
	# Inputs
	device = get_dict_from_yaml(device_yaml)  # yaml config file stored as dictionary
	layers = get_layers_from_yaml(device) 	  # a list of layer objects
	field_amp = device['wave']      		  # Electric field amplitude
	theta_i = field_amp['theta_i']  		  # Initial incident wave angle
	theta_f = field_amp['theta_f'] 			  # Final incident wave angle
	num_angles = field_amp['num_angles']  	  # Number of angles to sweep through
	angles = np.linspace(theta_i, theta_f, num_angles)
	
	# Initialize Wave class	
	min_wavelength = float(device['min_wavelength'])
	max_wavelength = float(device['max_wavelength'])
	num_points = int(device['num_points'])
	wave = Wave(min_wavelength, max_wavelength, num_points)
	wave.make_wavelengths()  # still in units of um

	# Interpolating downloaded index data so number of data points match.
	for layer in layers:
		layer.make_new_data_points(wave.wavelengths)

	# Make folder for simulation results
	device_name = device_yaml.split('/')[-1]  # Get filename without path or '.yaml'
	device_name = device_name[ : device_name.find('.yaml')]
	sim_foldername = device_name + '_' + str(min_wavelength) + '-' + str(max_wavelength) + 'um' + '_' + str(theta_i) + '-' + str(theta_f) + 'deg'
	sim_path = os.path.join(output_dir, sim_foldername)
	if not os.path.exists(sim_path):
		os.makedirs(sim_path)

	print("")
	num_cores = multiprocessing.cpu_count()
	print("CPU Core Count:", num_cores)
	pool = multiprocessing.Pool(num_cores)
	
	n = len(angles)
	pbar = tqdm(total=n)
	
	res = [pool.apply_async(perform_transfer_matrix, 
							args=(sim_path, angle, wave.wavelengths, layers, wave_type),
							callback = lambda _: pbar.update(1)) for angle in angles]
	results = [p.get() for p in res]
	pool.close()
	pool.join()
	pbar.close()
	print("")
	print("Wrote results to {}".format(sim_path))


def get_console_handler():
	console_handler = logging.StreamHandler(sys.stdout)
	console_handler.setFormatter(FORMATTER)
	return console_handler

def get_logger(args, logger_name):
	"""Configure logging."""
	logger = logging.getLogger(logger_name)
# 	logger.setLevel(logging.DEBUG)
	if args.debug:
		logger.setLevel(logging.DEBUG)
	else:
		logger.setLevel(logging.INFO)
	logger.addHandler(get_console_handler())
	logger.propagate = False
	return logger

def parse_arguments():
	parser = argparse.ArgumentParser()
	device_help = "Path for a yaml file from config_files describing a device."
	output_help = "Directory for transfer matrix results."
	pwave_help = "Boolean. Incident wave is p-wave."
	swave_help = "Boolean. Incident s-wave."
	units_help = "Choose the units for the output electric field. Default is micrometers."

	parser.add_argument('--debug', action='store_true', help="Enable debugging.")
	parser.add_argument("device", help=device_help)
	parser.add_argument("output", help=output_help)
	parser.add_argument('-p', '--pwave', help=pwave_help, action='store_true')
	parser.add_argument('-s', '--swave', help=swave_help, action='store_true')

	return parser.parse_args()


def main(args):

	logger.debug("Debugging enabled")
	logger.info("Start simulation")
	logger.info("Loading device parameters from {}".format(args.device))

	wave_type = 'mixed'
	if args.pwave:
		wave_type = 'p-wave'
		logger.info("Incident wave is p-wave.")
	elif args.swave:
		wave_type = 's-wave'
		logger.info("Incident wave is s-wave.")
	else:
		logger.info("Incident wave is mixed.")

	start_time = time.time()
	angle_resolved_multiprocess(args.device, args.output, wave_type)
	end_time = time.time()
	elapsed_time = np.round(end_time - start_time, 4)
	logger.info('Elapsed time: {} seconds'.format(elapsed_time))


if __name__ == '__main__':

	args = parse_arguments()
	logger = get_logger(args, __name__)
	main(args)
