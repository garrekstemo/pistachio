"""
Name: Transfer Matrix
Author: Garrek Stemo
Contact: garrekstemo@icloud.com
Affiliation: Nara Institute of Science and Technology
"""

import argparse
import csv
import importlib.resources as pkg_resources
import logging
import multiprocessing
import numpy as np
import os
import sys
from ruamel.yaml import YAML
import scipy.constants as const
from scipy import interpolate
import time
from tqdm import tqdm

from pistachio.data import refractive_index_data

yaml = YAML()
# FORMATTER = logging.Formatter("%(asctime)s — %(name)s — %(levelname)s - %(message)s")
FORMATTER = logging.Formatter("%(message)s")
LOG_FILE = 'transfer_matrix.log'


class Layer:
	"""
	Contains information and functions associated with a single, planar layer.

	Attributes
	----------
	material : string
		Name of the material comprising this layer.
	thickness : float
		Thickness of the layer (in meters).
	refractive_index : [...,...,...] 1d array
		Real part of the complex refractive index.
	extinction_coeff : [..., ..., ...] 1d array
		Imaginary part of the complex refractive index.
	"""
	def __init__(self, material: str='Air', thickness: float=1.0e-3, num_points: int=100):

		self.material = material
		self.thickness = thickness
		self.wavelengths = []
		self.refractive_index = []
		self.extinction_coeff = []
		self.set_complex_refractive(self.refractive_index, self.extinction_coeff)
		self.kx = []
		self.kz = []
		self.transfer_matrices = []

	def set_complex_refractive(self, n_real: list[float] = None, n_imag: list[float] = None) -> None:
		"""
		Sets the complex refractive index.

		Parameters
		----------
		n_real : [...,...,...] 1d array of floats
			List of real part of the refractive index for different wavelengths.
		n_imag : [...,...,...] 1 array of floats
			List of imaginary part of the refractive index for different wavelengths.
		
		Returns
		-------
		None
		"""
		n_complex = []
		for ii in range(len(n_real)):
			n_complex.append(n_real[ii] + 1j * n_imag[ii])
		self.complex_refractive = np.array(n_complex)

	def get_data_from_csv(self, refractive_filename: str) -> None:
		"""
		Extracts refractive index data associated with a layer material 
		from file downloaded from refractiveindex.info and sets wavelength,
		real refractive index, and extinction coefficient (if present).

		Parameters
		----------
		refractive_filename : str
			Name of the .csv file in data/refractive_index_data
			from which to extract refractive index data.

		Returns
		-------
		None

		Notes
		----
		The data from refractiveindex.info uses micrometers for wavelength units.
		This function assumes the file is left as is and converts to meters.
		"""
		with pkg_resources.path(refractive_index_data, refractive_filename) as params:
			params = os.path.abspath(params)

			wavelen = []
			n_real = []
			n_imag = []

			with open(params, 'r', encoding='utf-8') as csv_file:
				csvreader = csv.reader(csv_file)
				next(csvreader, None)
				for row in csvreader:
					wavelen.append(float(row[0]) * 10**-6)
					n_real.append(float(row[1]))
					try:
						n_imag.append(float(row[2]))
					except IndexError:
						n_imag.append(0.0)

			self.wavelengths = wavelen
			self.refractive_index = n_real
			self.extinction_coeff = n_imag
			self.set_complex_refractive(n_real, n_imag)

	def make_datapoints(self, wavelengths: list) -> None:
		"""
		Makes a new set of data points from user-defined num_points and max/min wavelength,
		and uses SciPy interpolation to match spacing and number of data points
		for refractive index so that array lengths are consistent between layers.

		SciPy's interpolate.interp1d uses x and y values to
		generate a function whose argument uses interpolation to find the value of new points.

			f(x) = y

		See the documentation for details.
		https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html#scipy.interpolate.interp1d

		Parameters
		----------
		wavelengths: [..., ..., ...] 1d array
			List of user-defined wavelengths of light (in meters) 
			via max, min wavelength in yaml config file.
	
		Returns
		-------
		None

		Notes
		-----
		Interpolation of the refractive index passes unit testing for values in the middle of the new 
		values for x (wavelengths) passed into the interpolation functions. Original refractive index
		and extinction coefficients are usually given to about three decimal places. Interpolated
		values sometimes preserve the original value to the full precision, but sometimes differ past 
		the first or second demical place. This is something the user ought to be aware of.
		The user should conduct independent tests of the accuracy of their own data.

		Furthermore, data near the edges of x are unstable and deviate a few percent from the reference data.
		When constructing new wavelength data using this function, DO NOT 
		use wavelengths outside of any reference refractive index data that you provide.
		The function will extrapolate that data instead of raising an error, and the
		new refractive indices may be unphysical.

		"""
		num_points = len(wavelengths)
		if len(self.refractive_index) == 1:
			self.refractive_index = np.full(num_points, self.refractive_index)
			self.extinction_coeff = np.full(num_points, self.extinction_coeff)
			self.set_complex_refractive(self.refractive_index, self.extinction_coeff)
			self.wavelengths = wavelengths
		else:
			f_n = interpolate.interp1d(self.wavelengths, self.refractive_index, fill_value='extrapolate')
			f_K = interpolate.interp1d(self.wavelengths, self.extinction_coeff, fill_value='extrapolate')

			self.refractive_index = f_n(wavelengths)
			self.extinction_coeff = f_K(wavelengths)
			self.set_complex_refractive(self.refractive_index, self.extinction_coeff)
			self.wavelengths = wavelengths

	def calculate_wavenumbers(self, n_complex: list, theta: float=0.0) -> None:
		"""
		Calculate wavenumbers for the dielectric for the given
		angular frequency and angle.

		Parameters
		----------
		n_complex : ndarray
			Array of complex refractive indices for a set of wavelengths.
		
		Returns
		-------
		None
		"""
		omega = 2 * np.pi * const.c / self.wavelengths
		self.kx = n_complex * omega / const.c * np.cos(theta)
		self.kz = n_complex * omega / const.c * np.sin(theta)

	def dynamical_matrix(self, n_complex: complex, theta: float, wave_type: str = 's-wave') -> list:
		"""
		Calculates the dynamical matrix for a layer.

		Parameters
		----------
		n_complex : complex
			Complex refractive index
		theta : float
			Incident angle of light with respect to normal in radians.
		wave_type : str
			Either s-wave or p-wave

		Returns
		-------
		D : numpy ndarray
			2x2 dynamical matrix.
		"""
		if wave_type == 's-wave':
			return np.array([[1, 1], [n_complex * np.cos(theta),-n_complex * np.cos(theta)]])
		
		elif wave_type == 'p-wave':
			return np.array([[np.cos(theta), np.cos(theta)], [n_complex, -n_complex]])

		else:
			print("Must specify 's-wave' or 'p-wave'.")

	def propagation_matrix(self, kx: float, d: float) -> np.ndarray:
		"""
		Calculates the propagation matrix through a layer
		(accumulation of phase for a plane wave propagating through homogeneous medium).

		Parameters
		----------
		kx : float
			Wavenumber associated with the incident wave.
		d : float
			Layer thickness in meters.

		Returns
		-------
		P : numpy ndarray
			Propagation matrix for the layer.
		"""
		return np.array([[np.exp(1j * kx * d), 0.0], [0.0, np.exp(-1j * kx * d)]])

	def calculate_transfer_matrix(self, n_complex: complex, lmbda: float, theta: float, wave_type: str) -> list:
		""" 
		Performs transfer matrix for this layer for a particular
		frequency and incidence angle of light.

		M(omega) = D_i * P_i * D_i^-1

		Parameters
		----------
		n_complex : complex number
			Complex refractive index for the given material and input wave.
		lmbda : float
			Wavelength of incoming light.
		theta : float
			Angle of incident light.
		wave_type : str
			Polarization, either 's-wave' or 'p-wave'.

		Returns
		-------
		matrices : list
			List containing dynamical matrix, propagation matrix,
			inverse dynamical matrix, and layer transfer matrix.

		Notes
		-----
		This really has to independent of the intrinsic attributes
		of the layer class, since we want to sweep through both angles and wavelengths
		of light. Incorporate this function into a larger loop when this sort of 
		sweep is desired. If you just want to calculate the transfer matrix
		for a single angle and wavelength, this is then still possible.

		"""
		omega = 2 * np.pi * const.c / lmbda
		kx = (n_complex * omega / const.c) * np.cos(theta)

		D_i = self.dynamical_matrix(n_complex, theta, wave_type)
		D_i_inv = np.linalg.inv(D_i)
		P_i = self.propagation_matrix(kx, self.thickness)
		M_i = np.linalg.multi_dot([D_i, P_i, D_i_inv])

		return [D_i, P_i, D_i_inv, M_i]

	def calculate_all_transfer_matrices(self, theta: float, wave_type: str) -> None:
		"""
		Calculate transfer matrices for this layer for every stored wavelength.

		Parameters
		----------
		theta : float
			Incident angle.
		wave_type : str
			's-wave' or 'p-wave' incident light polarization.

		Returns
		-------
		None

		Notes
		-----
		This is useful when you are ready to perform the total transfer matrix on
		the entire structure and number of wavelengths for each layer is the same.

		"""
		transfer_matrices = []
		for i, lmbda in enumerate(self.wavelengths):
			D_i, P_i, D_i_inv, M_i = self.calculate_transfer_matrix(self.complex_refractive[i], lmbda, theta, wave_type)
			transfer_matrices.append([D_i, P_i, D_i_inv, M_i])

		self.transfer_matrices = np.array(transfer_matrices)


class Structure:
	"""
	Structure class. A multilayered optical structure consisting of planar layer classes.

	Attributes
	----------
	layers : [..., ..., ...] list of layer objects
		List of all layers in the structure.

	"""
	def __init__(self):
		self.layers = []
		self.theta = []
		self.wavelengths = []
		self.wavevectors = []
		self.transfer_matrices = []

	def add_layer(self, layer: object) -> None:
		"""
		Append a new layer to the end of the list of layers.

		Parameters
		----------
		layer : layer object
			An instance of the Layer class.

		Returns
		-------
		None
		"""
		self.layers.append(layer)

	def delete_layer(self, layer_index: int) -> None:
		"""
		Remove a layer from the structure according to the given index.

		Parameters
		----------
		layer_index : int
			The index of an existing layer.
		
		Returns
		-------
		None
		"""
		self.layers.pop(layer_index)

	def load_yaml_config(self, config_file: str) -> dict:
		"""
		Load a structure from a yaml configuration file.
		See the provided templates for how this file should be structured.

		Parameters
		----------
		config_file : str
			The path to the yaml configuration file.
		
		Returns
		-------
		structure : dict
			A Python dictionary containing information about the structure,
			including the parameters the incident light waves,
			structure angle, and parameters for each layer.
 
		"""
		structure = dict()
		with open(config_file, 'r') as yml:
			structure = yaml.load(yml)
		return structure

	def load_struct_from_config(self, config_file: str) -> None:
		"""
		Load a structure stored in a .yaml config file and assigns layers
		as attributes of the structure.

		Parameters
		----------
		config_file : str
			Path to .yaml file
		Returns
		-------
		None

		"""
		struct_dict = self.load_yaml_config(config_file)

		num_layers = len(struct_dict['layers'])
		min_wl = float(struct_dict['min_wavelength'])
		max_wl = float(struct_dict['max_wavelength'])
		num_wl = int(struct_dict['num_points'])
		theta_i = float(struct_dict['wave']['theta_i'])
		theta_f = float(struct_dict['wave']['theta_f'])
		num_angles = int(struct_dict['wave']['num_angles'])
		polarization = str(struct_dict['wave']['polarization'])

		layers = []
		for i in range(num_layers):

			layer_str = 'layer' + str(i)
			layer = struct_dict['layers'][layer_str]
			layer_class = Layer()
			layer_class.material = str(layer['material'])
			layer_class.thickness = float(layer['thickness'])

			if "refractive_filename" in layer:
				params = layer['refractive_filename']
				if 'csv' in params:
					layer_class.get_data_from_csv(params)
				else:
					print("Cannot get refractive index data. Check file type (must be .csv).")

			elif "refractive_index" in layer:
				layer_class.refractive_index = np.array([layer['refractive_index']])
				layer_class.extinction_coeff = np.array([layer['extinction_coeff']])
				layer_class.set_complex_refractive(layer_class.refractive_index, layer_class.extinction_coeff)
				layer_class.wavelengths = np.array([layer['wavelength']])

			else:
				print("ERROR: Incorrect yaml config format. Reference default template.")

			layers.append(layer_class)

		self.layers = layers
		print('Initializing structure...')
		self.initialize_struct(theta_i, theta_f, num_angles, min_wl, max_wl, num_wl, polarization)

	def initialize_struct(self, theta_i: float, theta_f: float, num_angles: int, min_wl: float, max_wl: float, num_wl: int, wave_type: str) -> None:
		"""
		Initialize a multi-layered structure, including initial calculations
		and setting the number of data points for each layer.
		Initialize with first incident angle, theta_i.

		Use this after adding all desired layers.

		Parameters
		----------
		min_wl : float
			Starting (smallest) wavelength (in meters).
		max_wl : float
			Ending (largest) wavelength (in meters).
		num_points : int
			Number of data points to generate between max and min wavelengths

		Returns
		-------
		None

		"""
		self.wavelengths = np.linspace(min_wl, max_wl, num_wl)
		self.theta = np.deg2rad(np.linspace(theta_i, theta_f, num_angles))

		for layer in self.layers:
			layer.make_datapoints(self.wavelengths)
			layer.calculate_wavenumbers(layer.complex_refractive, theta_i)
			layer.calculate_all_transfer_matrices(theta_i, wave_type)

	def calculate_t_r(self, M: list):
		"""
		Calculate transmittance and reflectance coefficients and
		transmission and reflectance for a single layer.

		Parameters
		----------
		M : ndarray
			Transfer matrix. This can be for the whole structure, or just one or a few layers.
		
		Returns
		-------
		T : float
			Transmission
		R : float
			Reflection

		"""
		t = 1 / M[0,0]
		r = M[1,0] / M[0,0]
		
		R = np.real(np.abs(r)**2)
		T = np.real(np.linalg.det(M) * np.abs(t)**2)

		return T, R

	def calculate_all_t_r(self, transfer_matrices: list) -> list[float]:
		"""
		Calculate reflectance and transmittance for all stored transfer matrices,
		essentially constructing a transmittance and reflectance spectrum.

		Parameters
		----------
		None

		Returns
		-------
		T_spectrum : 1d numpy array
			Array of transmittance values at each stored wavelength.
		R_spectrum : 1d numpy array
			Array of reflectance values at each stored wavelength.
		"""
		T_spectrum = []
		R_spectrum = []

		for m in transfer_matrices:
			T_i, R_i = self.calculate_t_r(m)
			T_spectrum.append(T_i)
			R_spectrum.append(R_i)

		return np.array(T_spectrum), np.array(R_spectrum)

	def calculate_layer_matrices(self, theta: float, wave_type: str) -> None:
		"""
		Calculates the transfer matrices for each layer in the structure
		for the available wavelengths.

		Parameters
		----------
		theta : float
			Incident light angle.
		wave_type : str
			Light polarization, 's-wave' or 'p-wave'.

		Returns
		-------
		None

		"""
		for l in self.layers:
			l.calculate_all_transfer_matrices(theta, wave_type)

	def calculate_all_transfer_matrices(self, theta: float, wave_type: str) -> list:
		"""
		Calculate transfer matrices for every stored wavelength
		for the entire structure. Make sure the data are treated
		such that each layer has the same number of points. Calculated 
		transfer matrices are stored in a list as an attribute of the 
		Structure object, which can be used to calculate transmittance
		and reflectance spectra.

		Parameters
		----------
		theta : float
			Angle of incident light.
		wave_type : str
			Light polarization, either 's-wave' or 'p-wave'.
		
		Returns
		-------
		None
		"""
		self.calculate_layer_matrices(theta, wave_type)

		transfer_matrices = []
		for i in range(len(self.wavelengths)):

			M_build = np.identity(2, dtype=complex)
			D_0_inv = self.layers[0].transfer_matrices[i][2]
			D_f = self.layers[-1].transfer_matrices[i][0]

			for j in range(len(self.layers))[-2:0:-1]:
				M_j = self.layers[j].transfer_matrices[i][-1]
				M_build = np.matmul(M_j, M_build)

			M = np.linalg.multi_dot([D_0_inv, M_build, D_f])

			transfer_matrices.append(M)

		return transfer_matrices

	def one_job(self, theta: float, wave_type: str) -> list:
		"""
		Calculates transfer matrix, transmittance, and reflectance
		for a single angle. This is part of the angle-resolved calculation
		using Python's multiprocess module.

		Parameters
		----------
		theta : float
			Incidence angle.
		wave_type : str
			Polarization, either 's-wave' or 'p-wave'.
		
		Returns
		-------
		T_angle_resolved : list of list of floats
			List (for each angle) of list (for each wavelength) 
			of transmittance values.
		"""
		M_i = self.calculate_all_transfer_matrices(theta, wave_type)
		T_i, R_i = self.calculate_all_t_r(M_i)
		return T_i, R_i

	def angle_resolved_spectra(self, wave_type: str) -> list:
		"""
		Calculates transfer matrix for each provided incidence angle.
		Uses Python's multiprocess module to speed up calculation (useful
		for large number of angles).

		Parameters
		----------
		wave_type : str
			Polarization, either 's-wave' or 'p-wave'.

		Returns
		-------
		M_angle_resolved : list of lists
			List of list of transmittance and reflectance values
			for each light incident angle.
		"""
		num_cores = multiprocessing.cpu_count()
		pool = multiprocessing.Pool(num_cores)
		n = len(self.theta)
		pbar= tqdm(total=n)

		res = [pool.apply_async(self.one_job,
								 args=(angle, wave_type),
								 callback=lambda _: pbar.update(1)) for angle in self.theta]

		results = [p.get() for p in res]
		pool.close()
		pool.join()
		pbar.close()

		return results

	def print_structure(self) -> None:
		"""
		Prints the current structure configuration.

		Parameters
		----------
		None

		Returns
		-------
		None
		"""
		print('Structure Configuration')
		for i, layer in enumerate(self.layers):

			print("Layer {} :".format(i), layer.material,
				  "| d = " + str(int(layer.thickness * 10**9)) + " nm"
			)
	def make_sim_dir(self, out_path: str, sim_name: str=None) -> None:
		"""
		Make the name of the output directory for 
		angle-resolved transfer matrix simulation.
		Uses structure attributes if dir_name is None.
		The transmittance and reflectance data for each angle
		will be written to a separate file within this folder.

		Parameters
		----------
		out_path : str
			Path at which simulation directory will be created.
		sim_name : str
			Name of simulation directory.

		Returns
		-------
		None

		"""
		sim_path = None
		if sim_name == None:
			sim_name = ""
			for layer in self.layers:
				sim_name += layer.material + "_"
			theta_i = np.round(np.rad2deg(self.theta[0]))
			theta_f = np.round(np.rad2deg(self.theta[-1]))
			lmbda_i = self.wavelengths[0]
			lmbda_f = self.wavelengths[-1]
			param_str = str(theta_i) + '--' + str(theta_f) + 'deg_' + str(lmbda_i) + "--" + str(lmbda_f) + "lambda/"
			sim_name += param_str
			sim_path = os.path.join(out_path, sim_name)
		
		else:
			sim_path = os.path.join(out_path, sim_path)
		
		return sim_path


def write_tmm_results(theta: float, wavelengths: list[float], trans: list[float], refl: list[float], output_dir: str) -> None:
	"""
	Writes transmission, reflectance, results for all
	wavelengths specified in transfer matrix calculation to a csv file
	in the specified output directory.

	Parameters
	----------
	theta : float
		Angle (in degrees) associated with the incidence angle
		supplied for this calculation.
	wavelengths : list of floats
		List of wavelengths to write to file.
	trans : list of floats
		List of transmittance data.
	refl : list of floats
		List of reflectance data.
	output_dir : str
		Path to directory to save transfer matrix results.

	Returns
	-------
	None
	"""
	file_name = str(theta) + 'deg_.csv'
	output_file = os.path.join(output_dir, file_name)

	with open (output_file, 'w', encoding='utf-8', newline='') as out_file:

		filewriter = csv.writer(out_file, delimiter=',')
		header = ['Wavelength',
				  'Transmittance',
				  'Reflectance'
				  ]
		filewriter.writerow(header)
		i = 0
		while i < len(trans):
			row = [wavelengths[i], trans[i], refl[i]]
			filewriter.writerow(row)
			i+=1


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
	structure_help = "Path for a yaml file from config_files describing a multilayered structure."
	output_help = "Directory for transfer matrix results."
	pwave_help = "Boolean. Incident wave is p-wave."
	swave_help = "Boolean. Incident s-wave."

	parser.add_argument('--debug', action='store_true', help="Enable debugging.")
	parser.add_argument("structure", help=structure_help)
	parser.add_argument("output", help=output_help)
	parser.add_argument('-p', '--pwave', help=pwave_help, action='store_true')
	parser.add_argument('-s', '--swave', help=swave_help, action='store_true')

	return parser.parse_args()


def main(args):

	logger.debug("Debugging enabled")
	logger.info("Loading structure parameters from {}".format(args.structure))

	wave_type = 'None'
	if args.pwave:
		wave_type = 'p-wave'
		logger.info("Incident wave is p-wave.")
	elif args.swave:
		wave_type = 's-wave'
		logger.info("Incident wave is s-wave.")
	else:
		logger.info("Incident wave must be 'p-wave' or 's-wave'.")
		sys.exit()

	s = Structure()
	s.load_struct_from_config(args.structure)
	logger.info("Start simulation")
	start_time = time.time()
	results = s.angle_resolved_spectra(wave_type)
	end_time = time.time()
	elapsed_time = np.round(end_time - start_time, 4)
	logger.info('Calculation time: {} seconds'.format(elapsed_time))

	out_path = s.make_sim_dir(args.output)
	if not os.path.exists(out_path):
		os.mkdir(out_path)

	for i, theta in enumerate(results):
		t_results = theta[0]
		r_results = theta[1]
		write_tmm_results(s.theta[i], s.wavelengths, t_results, r_results, out_path)
	print("Wrote results to {}".format(out_path))


if __name__ == '__main__':

	args = parse_arguments()
	logger = get_logger(args, __name__)
	main(args)
