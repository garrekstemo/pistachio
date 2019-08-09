#! /anaconda3/bin/python

# Convention used
# Psi(x, t) = Psi_0 * exp(i(kx - wt))
# n = n' + i*n'', where n' is real part of refractive index and n'' is imaginary part.

import argparse
import sys
import csv
from itertools import tee
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy as sp
import scipy.constants as sc
import scipy.interpolate
from ruamel_yaml import YAML
import pdb

c = sc.c  # speed of light
h = sc.h  # planck's constant
yaml = YAML()


class Light:

	def __init__(self, wavelength):
		
		self.wavelength = wavelength
# 		self.omega = 2*np.pi*c / wavelength  # angular frequency
# 		self.freq = c / wavelength  # frequency
# 		self.k = 2*np.pi / wavelength  # wavenumber
# 		self.energy_J = h*c / wavelength  # energy in Joules
# 		self.energy_ev = self.energy_J / sc.eV  # energy in electron-volts
		self.amplitude = []
		self.trans_amplitude = np.array([1, 0])
		self.matrices = []
		
	def omega(self):
		return 2*np.pi*c / self.wavelength

	def energy_J(self):
		return h*c / self.wavelength
	
	def energy_ev(self):
		return h*c / self.wavelength / sc.eV
		
	def wavenumber(self):
		return 2*np.pi / self.wavelength
		
	def frequency(self):
		return c / self.wavelength


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
		self.waves = {}
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
				try:
					K = float(row[2])
					self.extinct.append(K)
				except IndexError:
					self.extinct.append(0.)
				
	def get_data_from_txt(self, path):
		"""Used for filmetrics.com data"""
		with open(path, 'r') as params:
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
					 
		# s-wave dynamical matrix
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
			print("ERROR: Incorrect yaml config format. Reference default template.")
			
		layers.append(layer_class)

	return layers


<<<<<<< Updated upstream
=======
def get_beam_profile(beam_csv):
	"""Jasco FTIR has weird non utf-8 characters in the footer and partway through data.
	This annoyance is kinda dealt with."""
	
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
					pass
		except StopIteration:
			pass
		
	return wavenum, field

def attenuation_coefficient(beam_intensity, x):
	"""Takes beam intensity as array"""
	
	dI = np.gradient(beam_intensity, x)
	
	return dI / beam_intensity
	
def kramers_kronig(alpha):
	"""Take attenuation coefficient. Returns real part of index of refraction."""
	
	


>>>>>>> Stashed changes
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

def build_matrix_list(wave_obj, theta, bound1, bound2, layers):
	"""TODO: Make this part of the Light object?
	Makes a wavelength object and bounds. Outputs a list of matrices that will
	make the transfer matrix for that wavelength of light propagating 
	through the device."""

	wave_key = wave_obj.wavelength   # Keep this key stored to retrive index data
	wave_obj.wavelength = wave_obj.wavelength * 10**-6  # used for computation
	omega = wave_obj.omega()
# 	print("build matrix, wave key", wave_key)
# 	print("build matrix, wave obj wl", wave_obj.wavelength)
	matrices = []
# 	lmbda = wave_obj[idx] * um
# 	light = Light(lmbda)

	n0 = 0
	if isinstance(bound1.index, np.ndarray):
		n0 = bound1.waves[wave_key]
	else:
		n0 = bound1.index + 1j*bound1.extinct
	D0 = bound1.dynamical_matrix(n0, theta)
	D0inv = np.linalg.inv(D0)
	matrices.append(D0inv)

	for layer in layers:
# 		n  = layer.index[idx] + 1j*layer.extinct[idx]
		n = layer.waves[wave_key]
		kx = layer.wavenumber(n, omega, theta)[0]

		D = layer.dynamical_matrix(n, theta)
		Dinv = np.linalg.inv(D)
		P = layer.propagation_matrix(kx)
		matrices.extend([D, P, Dinv])

	# Add transmission region medium to matrix list
	ns = 0
	if isinstance(bound2.index, np.ndarray):
		ns = bound2.waves[wave_key]
	else:
		ns = bound2.index + 1j*bound2.extinct
	Ds = bound2.dynamical_matrix(ns, theta)
	matrices.append(Ds)

	return matrices


def field_amp(matrix_arr, As_, Bs_):
	"""TODO: Make this part of the Light object?
	Takes transfer matrix, final magnitude of EM wave, As, Bs.
	Outputs E-field amplitude for given matrix. This approach starts calculations from
	the incident wave."""
	
	Elrev = []  #reversed order field amplitudes
	Es = np.array([As_, Bs_])
	Elrev.append(Es)
	i = 1
	W = 3
	L = len(matrix_arr)

	while W*i < L:
		slice = matrix_arr[L-W*i:]
		M = np.linalg.multi_dot(slice)
		amp = np.dot(M, Es)
		Elrev.append(amp)
		i+=1

	first_layer = 2
	M = np.linalg.multi_dot(matrix_arr[:first_layer])
	amp = np.dot(M, Es)
	Elrev.append(amp)

	# Put amplitudes in order from first to last layer.
	El = []
	for E in reversed(Elrev):  
		El.append(E)
	return El
	

def set_bound(device_, bound_name):
	"""Takes device dictionary and name of bounding material as a string"""
	
	bound_msg = "Boundary material name is not a string."
	assert isinstance(bound_name, str), bound_msg
	
	num_points = device_['num_points']
	MIN_WL = device_['min_wl']
	MAX_WL = device_['max_wl']
	name = device_[bound_name]['material']
	bound = Layer(name, num_points, MIN_WL, MAX_WL)
	param_path = 'param_path'
	
	if param_path in device_[bound_name]:
		params = device_[bound_name][param_path]
		bound.get_data_from_csv(params)
		bound.make_new_data_points()
	else:
		bound.index = device_[bound_name]['index']
		bound.extinct = device_[bound_name]['extinction']
# 		bound.make_new_data_points()
	
	return bound
		

# ===== Fresnel equation functions below not used so much for now ===== #

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

# def main_loop():
	


def main():
	# Inputs
	device = get_dict_from_yaml(args.device)  # yaml config file stored as dictionary
	layers = get_layers_from_yaml(device)  # a list of layer objects

	# If zero, p-wave and s-wave should yield same transmission
	em_wave = device['wave']
	inc_ang = em_wave['theta_in']
	
	# TODO: make code work for bounds with frequency dependent refractive indices
	bound1 = set_bound(device, 'bound1')
	bound2 = set_bound(device, 'bound2')
	
	print("{} layers".format(device['num_layers']))
	for layer in layers:
		# interpolating from downloaded index data, making new data points
		print(str(layer.material) + ", d=" + str(int(layer.thickness*10**9)) + "nm")
		layer.make_new_data_points()

	# We should separate this from the Layer class. 
	wavelens = layers[0].wavelength  # still in units of um
	all_of_the_lights = []
	
	# Outputs
	T = []
	R = []
	A = []
	E_amps = []
	intensity = []
# 	beams = "/Users/garrek/projects/raw_data/190731/0.csv"
# 	k, efield = get_beam_profile(beams)
	
# 	alpha = attenuation_coefficient(efield, k)

	for lm in wavelens:
		lmbda = Light(lm)
		all_of_the_lights.append(lmbda)
		M = build_matrix_list(lmbda, inc_ang, bound1, bound2, layers)
		lmbda.matrices = M
		TM = np.linalg.multi_dot(M)
		trns = transmittance(TM).real
		refl = reflectance(TM).real
		T.append(trns)
		R.append(refl)
		abso = 1 - trns - refl
		A.append(abso)

# 	for lm in all_of_the_lights:
# 		matrices = lm.matrices
# 		As, Bs = lm.trans_amplitude
# 		E_amps.append(field_amp(matrices, As, Bs))

# 	for amp in E_amps:
# 		E = amp[0][1]
# 		I = E*np.conj(E)
# 		intensity.append(I.real)


	# ===== Make Plots ===== #
	
	print("_"*50)
	print("Generating plots...")
	fig, axs = plt.subplots(3, 1, sharex=True)
	gs1 = gridspec.GridSpec(3, 1)
	gs1.update(wspace=0.025, hspace=0.005)

	ax = axs[0]
	ax.plot(wavelens, T, 'b-', label="simulation")
	if args.trns:
		ax.plot(wl_T_data, T_data, linestyle="dashed", color='#FF5733', label="downloaded data")
	ax.set_ylabel('Transmittance %', fontsize=16)
	ax.tick_params(axis='both', labelsize=16)
# 	ax.set_xlabel('wavelength ($\mu$m)', fontsize=20)
	ax.set_xlim(1, 10)
# 	ax.legend()
	
	
	ax = axs[1]
	ax.plot(wavelens, R, 'b-', label="calculated data")
	if args.refl:
		ax.plot(wl_R_data, R_data, linestyle="dashed", color='#FF5733', label="downloaded data")
	ax.set_ylabel('Reflectance %', fontsize=16)
	
	ax = axs[2]
	ax.plot(wavelens, A, 'b-', label="calculated data")
	if args.absorp:
		ax.plot(wl_A_data, A_data, linestyle="dashed", color="#FF5733", label="downloaded data")
	ax.set_ylabel('Absorptance %', fontsize=16)
	ax.set_xlabel('wavelength ($\mu$m)', fontsize=16)
	
	title = "Etalon with 10 nm Au film, 10 micron air gap"
# 	plt.suptitle(title, fontsize=24)
	plt.subplots_adjust(top=0.9)
# 	plt.tight_layout()
	plt.show()
	



# ===== Plotting FTIR data ===== #
# 	x_file = "/Users/garrek/Desktop/fpi0xwave.txt"
# 	y_file = "/Users/garrek/Desktop/trans0.txt"
# 	x_data = []
# 	y_data = []
# 	
# 	with open(x_file, 'r') as xf:
# 		lines = xf.readlines()
# 		for line in lines:
# 			x_data.append(float(line))
# 	with open(y_file, 'r') as yf:
# 		lines = yf.readlines()
# 		for line in lines:
# 			y_data.append(float(line))
# ============================== #

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	device_help = "path for a yaml file from config_files describing a device"
	reflectance_help = "path for reflectance data downloaded from filmetrics"
	transmittance_help = "path for transmittance data downloaded from filmetrics"
	absorptance_help = "path for absorptance data downloaded from filmetrics"
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
	parser.add_argument('-A', '--absorp', help=absorptance_help)
	debug_mode = parser.parse_args().debug
	if debug_mode:
		print("Debug mode\n")
	args = parser.parse_args()
	
	print("\nStart simulation\n")

	# Store a bunch of downloaded data
	if args.trns:
		wl_T_data, T_data = reference_data(args.trns)
	if args.refl:
		wl_R_data, R_data = reference_data(args.refl)
	if args.absorp:
		wl_A_data, A_data = reference_data(args.absorp)
	else:
		print("Using no reference data.\n")

	main()
