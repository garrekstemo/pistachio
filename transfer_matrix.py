#!/usr/bin/env python

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
# FORMATTER = logging.Formatter("%(asctime)s — %(name)s — %(levelname)s - %(message)s")
FORMATTER = logging.Formatter("%(message)s")
LOG_FILE = 'transfer_matrix.log'

class Wave:

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

	def __init__(self, material, num_points, min_wl=0.0, max_wl=0.0, thickness=0.0):
		self.material = material
		self.thickness = thickness
		self.num_points = num_points
		self.min_wl = min_wl  # starting wavelength
		self.max_wl = max_wl  # ending wavelength
		self.wavelengths = []  # array of free space wavelengths
		self.refractive_index = []  # array of refractive indices
		self.extinction_coeff = []  # array of extinction coefficients
		self.complex_refractive = {}  # Needs to be curly braces

	def __repr__(self):
		a = "{} \n".format(self.material)
		b = "thickness: {}\n".format(self.thickness)
		c = "wavelengths: {}\n".format(self.wavelengths)
		d = "refractive index: {}\n".format(self.refractive_index)
		e = "extinction coefficient: {}\n".format(self.extinction_coeff)
		return a+b+c+d+e

	def get_data_from_csv(self, index_path):
		"""Used for refractiveindex.info data. This site uses um for wavelength units."""
		with open(index_path, 'r') as params:
			reader = csv.reader(params)
			next(reader, None)
			for row in reader:
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

				self.wavelengths.append(wl)
				self.refractive_index.append(n)
				if line[2]:
					K = float(line[2])
					self.extinction_coeff.append(K)

# 	def set_wavelengths(self):
# 		"""Creates a new set of linearly-spaced wavelengths from start and
# 		   endpoint wavelengths and number of points specified in yaml config file."""
# 		new_wl = np.linspace(self.min_wl, self.max_wl, num=self.num_points, endpoint=True)
# 		return new_wl

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
	"""Takes device dictionary and outputs all layer objects as a list."""
	val1 = 'layers'
	val2 = 'num_points'
	val3 = 'min_wavelength'
	val4 = 'max_wavelength'

	num_layers = len(device_dict[val1])
	num_points = int(device_dict[val2])

	logger.info("Num points: {}".format(num_points))
	logger.info("Number of layers: {}".format(num_layers))

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
		logger.info(str(layer_class.material) + ", d=" + str(int(layer_class.thickness*10**9)) + "nm")

		if "param_path" in layer:
			params = layer['param_path']
			if 'txt' in params:
				layer_class.get_data_from_txt(params)
			elif 'csv' in params:
				layer_class.get_data_from_csv(params)
		elif "index" in layer:
			layer_class.refractive_index = layer['index']
			layer_class.extinction_coeff = layer['extinction']
			layer_class.wavelengths = layer['wavelength']
		else:
			logger.debug("ERROR: Incorrect yaml config format. Reference default template.")

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
		pol_msg = "\nNo wave acceptable polarization passed in.\nMixed-wave not yet supported.\nSee --help for more info. Exiting...."
		logger.debug(pol_msg)
		sys.exit()


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


def output_TRA(angle, output_dir, rows):
	"""Writes transmission, reflectance, abosorbance data to csv file"""

	wavelens, trans, refl, absor = rows
# 	output = os.path.abspath(args.output)
	file_name = 'deg' + str(angle)
	output_file = os.path.join(output_dir, file_name)
	with open (output_file, 'w') as out_file:
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
# 	logger.info("Wrote results to {}".format(output_file))
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
	logger.info("Wrote field output to {}".format(field_output))

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


def main_loop(device_yaml, output_dir, wave_type):
	"""
	Inputs: yaml file containing information about device and incident radiation.
	Outputs: Executes transfer matrix and other functions
	   		 Writes output file.
	"""

	# Inputs
	device = get_dict_from_yaml(device_yaml)  # yaml config file stored as dictionary
	layers = get_layers_from_yaml(device)  # a list of layer objects
	em_wave = device['wave']
	theta_i = em_wave['theta_i']  # Initial incident wave angle
	theta_f = em_wave['theta_f']  # Final incident wave angle
	num_angles = em_wave['num_angles']  # Number of angles to sweep through
	angles = np.linspace(theta_i, theta_f, num_angles)
	
	# Initialize Wave class
	
	min_wavelength = float(device['min_wavelength'])
	max_wavelength = float(device['max_wavelength'])
	num_points = int(device['num_points'])
	wave = Wave(min_wavelength, max_wavelength, num_points)
	wave.make_wavelengths()  # still in units of um

	logger.info('theta_i: {}, theta_f: {}, num angles: {}'.format(theta_i, theta_f, num_angles))
	for layer in layers:
		# interpolating from downloaded index data so number of data points match.
		#TODO: Am I doing this correctly?
		layer.make_new_data_points(wave.wavelengths)


	#TODO: Don't hard-code this directory when testing is done.
	sim_folder = os.path.join(output_dir, 'testing_sim_folder')
	if not os.path.exists(sim_folder):
		os.makedirs(sim_folder)

	for angle in angles:

		#TODO: Organize outputs better cuz there'll probably be a lot of them.
		# Outputs
		T = []
		R = []
		A = []
		E_amps = []
		intensity = []

	# 	alpha = attenuation_coefficient(efield, k)
		for lmbda in wave.wavelengths:

			M = build_matrix_list(lmbda, angle, layers, wave_type)
			TM = np.linalg.multi_dot(M)
			trns = transmittance(TM).real
			refl = reflectance(TM).real
	# 		field = field_amp(M, device["wave"]['A0'], device['wave']['B0'])
	# 		logging.info('Only using forward-propagating field values')
	# 		field = [f[0] for f in field]

	# 		E_amps.append(field)
			T.append(trns)
			R.append(refl)
# 			abso = 1 - trns - refl
			A.append(1 - trns - refl)

		#Write everything to a csv file
		output_TRA(angle, sim_folder, [wave.wavelengths, T, R, A])
	# 	output_field_profile(wavelens, layers, E_amps)


def get_console_handler():
	console_handler = logging.StreamHandler(sys.stdout)
	console_handler.setFormatter(FORMATTER)
	return console_handler

# def get_file_handler():
# 	file_handler = TimedRotatingFileHandler(LOG_FILE, when='midnight')
# 	file_handler.setFormatter(FORMATTER)
# 	return file_handler

def get_logger(args, logger_name):
	"""Configure logging."""
	logger = logging.getLogger(logger_name)
# 	logger.setLevel(logging.DEBUG)
	if args.debug:
		logger.setLevel(logging.DEBUG)
	else:
		logger.setLevel(logging.INFO)
	logger.addHandler(get_console_handler())
# 	logger.addHandler(get_file_handler())
	logger.propagate = False
	return logger

def parse_arguments():
	parser = argparse.ArgumentParser()
	device_help = "path for a yaml file from config_files describing a device"
	output_help = "path and name for output data file (must be .csv)"
	pwave_help = "boolean, calculate for p-wave"
	swave_help = "boolean, calculate for s-wave"

	parser.add_argument('--debug', action='store_true', help="Enable debugging.")
	parser.add_argument("device", help=device_help)
	parser.add_argument("output", help=output_help)
	parser.add_argument('-p', '--pwave', help=pwave_help, action='store_true')
	parser.add_argument('-s', '--swave', help=swave_help, action='store_true')
# 	parser.add_argument('-T', '--angles', action='store_true')

	return parser.parse_args()


def main(args):

	logger.debug("Debugging enabled")

	# NOTICE: USER NO LONGER HAS TO SPECIFY A FILE NAME

# 	output_message = "Output file must be .csv"
# 	assert args.output[-4:] == ".csv", output_message
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
	main_loop(args.device, args.output, wave_type)
	end_time = time.time()
	elapsed_time = np.round(end_time - start_time, 4)
	logger.info('elapsed time: {} seconds'.format(elapsed_time))


if __name__ == '__main__':

	args = parse_arguments()
	logger = get_logger(args, "TMM Logger")
	main(args)
