#! /anaconda3/bin/python


import argparse
import os
import re
import sys
import csv
import numpy as np
from pathlib import Path
from scipy import optimize
from scipy.interpolate import interp1d
from scipy import constants
from scipy import signal
from scipy import signal
import lmfit as lm
import matplotlib.pyplot as plt
from ruamel_yaml import YAML
import pdb
import convert_unit
import pmath

yaml = YAML()


# ========== Get paramaters and data from user inputs and files ========== #

def get_bounds_from_yaml(yaml_config):
	"""Takes a yaml config file and gets the upper and lower bound
		used for truncating and Lorentzian fitting later."""
	with open(yaml_config, 'r') as yml:
		config = yaml.load(yml)
		
	lower = config['bounds']['lower']
	upper = config['bounds']['upper']
	return [lower, upper]


def get_output_path_from_yaml(yaml_config):
	"""Gets data output path from yaml config file.
	   This saves the user a command-line argument.
	   You're welcome.'"""	
	   
	with open(yaml_config, 'r') as yml:
		config = yaml.load(yml)
	return config['data']['output']


def get_initial_from_yaml(yaml_config):
	"""Takes yaml config file and gets initial guesses
		for nonlinear least squares fit of dispersion curves later."""
		
	with open(yaml_config, 'r') as yml:
		config = yaml.load(yml)
	
	initial_units = config['least_squares_guesses']['units']
	E_0_guess = config['least_squares_guesses']['E_cav_0']
	Rabi_guess = config['least_squares_guesses']['Rabi_splitting']
	n_guess = config['least_squares_guesses']['refractive_index']
	E_vib_guess = config['least_squares_guesses']['E_exc']
	
	initial = [E_0_guess, E_vib_guess, Rabi_guess]

	if initial_units.lower() == 'ev':
		initial = convert_unit.set_units(initial, initial_units.lower(), 'wn')
	
	initial.append(n_guess)

	return initial, initial_units
	

def get_FTIR_data(spectral_data, data_format=None):
	"""
	#TODO: This has been renamed. Remember to adjust everything else downstream.
	Retrieve csv spectral data and store in array."""

	wavenumber_list = np.array([])
	intensity_list = np.array([])
	with open(spectral_data, 'r', errors='ignore') as spectrum:
		csvreader = csv.reader(spectrum)
		num_header_rows = 19   # Yeah this is hard coded. Not so great.
		for s in range(num_header_rows):
			next(csvreader, None)
		for row in csvreader:
			if row:
				wavenum = float(row[0])
				inten = float(row[1])
				wavenumber_list = np.append(wavenumber_list, wavenum)
				intensity_list = np.append(intensity_list, inten)
			else:
				break
	return wavenumber_list, intensity_list
	
def get_tmm_data(datafile):
	"""
	Retrieve data from transfer matrix simulation.
	"""


def get_angle_data_from_dir(directory, convert_units=None, data_format=None):
	"""
	Extracts angle-resolved and absorbance data from each file in the supplied directory.
	x-axis data is assumed to be in cm-1, but you can convert to other units here.
	
	convert_units: This should be a tuple containing the input units and the desired output units
				   Available units are in the pmath Convert class.
				   Example: convert_units = ('cm-1', 'um')
	"""

	abs_wavenum = np.array([]) 			 # absorbance wavenumbers
	abs_intensity = np.array([])			 # absorbance intensities
	angle_data = []  # List with structure [[angle, [wavenumbers, intensities]], [...], ...]

	#TODO: spectrum = spectrum.lower()
	deg_str = 'deg'
	abs_str = 'Abs'

	for spectrum in os.listdir(directory):
		if spectrum.endswith('.csv'):
			spec_file = str(directory) + '/' + str(spectrum)
			if deg_str in spectrum:
				# Get the angle of the measurement from file string
				deg_s = spectrum.find(deg_str) + len(deg_str)
				deg_e = spectrum.find('_', deg_s)
				deg = float(spectrum[deg_s:deg_e])
				x_data, intensity = get_FTIR_data(spec_file)

				if convert_units:
					x_data = pmath.set_units(x_data, convert_units[0], convert_units[1])

				angle_data.append([deg, x_data, intensity])

			# Get absorbance data if the file exists (it should always exist)
			if abs_str in spectrum:
				abs_wavenum, abs_intensity = get_FTIR_data(spec_file)

	absorbance_data = np.stack((abs_wavenum, abs_intensity))
	angle_data.sort()
	
	return angle_data, absorbance_data


def get_param_from_string(string, separator):
	"""Use regular expressions to find character in string
	   with repeated separation character (like '/' or '_')"""
	   
	str_loc = [m.start() for m in re.finditer(separator, string)]
	params = []
	p = 0
	num_locs = len(str_loc)
	len_sep = len(separator)

	count = 0
	while count < num_locs:

		if count == 0:
			first = 0
			second = str_loc[count]
			third = str_loc[count+1]
			p = string[first:second]
			params.append(p)
			p = string[second:third]
			
		elif str_loc[count] == str_loc[-1]:
			first = str_loc[count]
			p = string[first:]

	return angle_data, absorbance_data


def get_sample_params(directory):
	"""Gets parameters for each angle from angle-resolved directory and file name.
	   Returns concentration (M = mol/liter), solute name, solvent name."""

	sample_name = directory.split('/') #last element is angle directory
	sample_name = list(filter(None, sample_name))[-1]
	params = sample_name.split('_')

	if 'in' in params:
		params.remove('in')
	return sample_name, params


def truncate(xdata, ydata, bound1, bound2):
	"""Truncate data to isolate desired peaks for fitting"""
	i = 0
	MAX = len(xdata)
	zipped = zip(xdata, ydata)
	sorted_data = sorted(zipped, key = lambda t: t[0])
	xdata, ydata = zip(*sorted_data)
	xdata = np.asarray(xdata)
	ydata = np.asarray(ydata)
	lower_pt = 0
	upper_pt = 0
	
	for i, val in enumerate(xdata):
		if np.round(val, 0) == float(bound1):
			lower_pt = i

		if np.round(val, 0) == float(bound2):
			upper_pt = i
	
	return xdata[lower_pt:upper_pt], ydata[lower_pt:upper_pt]


# ============== Write results to files ============== #

def write_angle_spec_to_file(angle_data_list, sample_name, out_path):
	"""Takes angles list, wavenumber list, and intensity list.
		Writes angle-resolved data to csv file with
		wavenumber in leftmost column (x-axis), and intensities
		for each degree (y-axes)."""

	angles = []
	wavenumbers = angle_data_list[0][1]
	num_wavenums = len(wavenumbers)

	angle_res_file = sample_name + '_angle-resolved_spectra.csv'
	output = os.path.join(os.path.abspath(out_path), angle_res_file)
	with open(output, 'w') as out_file:
		filewriter= csv.writer(out_file, delimiter=',')
		header = ['Wavenumber (cm-1)']
		for a in angle_data_list:
			angle = a[0]
			header.append('Int deg' + str(angle) + ' (arb)')
		filewriter.writerow(header)

		i = 0
		# BIG ASSUMPTION: all data sets have same number of wavenumbers with the same values
		while i < num_wavenums:
			# Cycle through data points, first set up wavenumbers
			row = [wavenumbers[i]]

			for item in angle_data_list:
				# Cycle through data point for each degree, make a row
				row.append(item[2][i])

			filewriter.writerow(row)
			i+=1

	print('Wrote angle-resolved spectra results to {}\n'.format(output))


def write_dispersion_to_file(angles, E_up, E_lp, E_vib, E_cav, sample_name, out_path):
	"""Takes angles and energies (wavenumbers) for upper and lower
	   polaritons, vibrational energy, calculated cavity mode
	   and writes all that data to a csv file.
	   REMEMBER: change plot_dispersion in plots.py if output file format changes."""

	dispersion_file = sample_name + '_dispersion.csv'
	output = os.path.join(os.path.abspath(out_path), dispersion_file)
	with open(output, 'w') as f:
		filewriter = csv.writer(f, delimiter=',')
		header = ['Angle (deg)',
				  'UP Wavenumber (cm-1)', 'LP Wavenumber (cm-1)',
				  'Vibration mode (cm-1)', 'Cavity mode (cm-1)']
		filewriter.writerow(header)

		i = 0
		while i < len(angles):
			row = [angles[i], E_up[i], E_lp[i], E_vib, E_cav[i]]
			filewriter.writerow(row)
			i+=1

	print('Wrote dispersion results to {}'.format(output))


def write_peakfit_residuals_to_file(wavenum, residuals, angle, outpath):
	"""Writes dispersion residuals to file to evaluate the goodness of fit
	   for peak fitting algorithm."""
	
	# Make a directory for the files we are about to make
	residual_dir = os.path.join(outpath, 'residuals')
	Path(residual_dir).mkdir(parents=True, exist_ok=True)
	residual_filename = 'residuals_{}deg.csv'.format(angle)
	output = os.path.join(residual_dir, residual_filename)

	with open(output, 'w') as rf:
		filewriter = csv.writer(rf, delimiter=',')
		i = 0
		while i < len(wavenum):
			row = [wavenum[i], residuals[i]]
			filewriter.writerow(row)
			i+=1


def write_splitting_fit_to_file(least_squares_results, sample_name, units, out_path):
	"""Takes results of nonlinear least squares fit for Rabi splitting
	   (E_cav_0, Rabi, n, E_vib) and writes them to a file for 
	   later plotting and whatnot."""
	
	E_cav_0, E_vib, Rabi, n_ = least_squares_results
	least_squares_file = sample_name + '_splitting_fit.csv'
	output = os.path.join(os.path.abspath(out_path), least_squares_file)
	
	with open(output, 'w') as f:
		filewriter = csv.writer(f, delimiter=',')
		filewriter.writerow(['units', units])  #TODO: What units to write?
		filewriter.writerow(['E_cav_0', E_cav_0])
		filewriter.writerow(['Rabi', Rabi])
		filewriter.writerow(['n', n_])
		filewriter.writerow(['E_vib', E_vib])

	print('Wrote results of nonlinear least squares to {}'.format(output))


def error_f(x, theta, Elp_data, Eup_data):
	"""Here, we write the eigenvalues of the energy Hamiltonian for 
	   coupled cavity-vibration system. We minimize this with experimental data."""

	En = pmath.coupled_energies(theta, *x, branch=0)
	Ep = pmath.coupled_energies(theta, *x, branch=1)
	
	err1 = En - Elp_data
	err2 = Ep - Eup_data	
	err = np.concatenate((err1, err2))
			
	return err


def optimize_df(x, theta, E_up, E_lp):
	"""The derivative w.r.t. each variable for optimize_f function to
	   compute Jacobian in nonlinear least squares fit."""
	
	n = len(2*theta)
	m = len(x)
	jacobian = np.zeros((n, m))   
	
	E_0, E_vib, Rabi, n_eff = x
	E_cav = E_0 / np.sqrt(1 - np.sin(theta)**2 /n_eff**2)  #Can probably use the function
	
	dEc_dE0 = 1 / np.sqrt(1 - (np.sin(theta)**2 / n_eff**2))
	dEc_dn = E_0 * np.sin(theta)**2 / (n_eff*n_eff*n_eff * np.power((1 - (np.sin(theta)**2 / n_eff**2)), 1.5))
	root = np.sqrt(4*Rabi**2 + (E_vib - E_cav)**2)

	# Negative solution	
	dE_vib_neg = 0.5 - (0.5*(E_vib - E_cav) / root)
	dRabi_neg = -2*Rabi / root
	dE_0_neg = (0.5 + 0.5*(E_vib - E_cav) / root) * dEc_dE0
	dn_neg = (-0.5 - 0.5*(E_vib - E_cav) / root) * dEc_dn

	# Positive solution
	dE_vib_pos = 0.5 + (0.5*(E_vib - E_cav) / root)
	dRabi_pos = 2*Rabi / root
	dE_0_pos = (0.5 - 0.5*(E_vib - E_cav) / root) * dEc_dE0
	dn_pos = (-0.5 + 0.5*(E_vib - E_cav) / root) * dEc_dn

	dE_vib = np.concatenate((dE_vib_neg, dE_vib_pos))
	dRabi = np.concatenate((dRabi_neg, dRabi_pos))
	dE_0 = np.concatenate((dE_0_neg, dE_0_pos))
	dn_eff = np.concatenate((dn_neg, dn_pos))
	
	jacobian[:,0] = dE_0
	jacobian[:,1] = dRabi
	jacobian[:,2] = dn_eff
	jacobian[:,3] = dE_vib

	return jacobian


def splitting_least_squares(initial, theta, Elp, Eup):
	"""Takes initial guesses: [E_cav_0, E_vib, Rabi, n].
	   Takes angles (in degrees), and experimental upper and lower polariton data.
	   Returns nonlinear least squares fit."""

	theta_rad = [np.pi/180*a for a in theta]
	optim = optimize.least_squares(error_f, jac='3-point',
								   x0=initial,
								   args=(theta_rad, Elp, Eup))
	return optim


def parse_args():
	parser = argparse.ArgumentParser()

	spectrum_help = "single csv file or directory containing spectral data."
	fit_function_help = "String that determines which fitting function to use when doing peak \
						 fitting of polariton spectral data. \
						 Options include 'single_peak' and 'double_peak' \
						 for two single-peak Lorentzian fuctions and one double-peak \
						 Lorentzian function, respectively."
	config_help = "Yaml file. Truncates data to include only polariton peaks; \
					sets least squares fit initial guesses for: \
					0 degree incidence cavity mode energy, \
					Rabi splitting value, \
					vibrational excitation energy; sets output file path. \
					See ReadMe for more information."

	pol_help = "Boolean. Lorentzian fit for a single spectrum file."
	abs_help = "Boolean. Lorentzian fit for absorbance single peak data."
	angle_help = "Boolean. Directory contains angle-tuned data."


	parser.add_argument('spectral_data', help=spectrum_help)
	parser.add_argument('--config', help=config_help)
	parser.add_argument('-P', '--polariton', action='store_true', help=pol_help)
	parser.add_argument('-A', '--absorbance', action='store_true', help=abs_help)
	parser.add_argument('-T', '--angle', action='store_true', help=angle_help)
	parser.add_argument('-C', '--concentration', action='store_true')
	parser.add_argument('fit_function', type=str, nargs='*', help=fit_function_help)

	return parser.parse_args()


def main():
	args = parse_args()
	spectral_data = args.spectral_data

	if args.config:
	
		config_params = args.config
		bounds = get_bounds_from_yaml(config_params)
		bounds.sort()
		output_path = get_output_path_from_yaml(config_params)

	if args.fit_function:
		fit_func = args.fit_function[0]

	if args.polariton:
		print('Fitting double-peak Lorentzian')
		x, y, fit = polariton_fitting(spectral_data)
		plot_polariton(x, y, fit)

	elif args.absorbance:
		print('Fitting single-peak Lorentzian')
		k, I = get_data(spectral_data)
		lor = absorbance_fitting(k, I)
		plot_absorbance(k, I, lor)

	elif args.angle:

		print('Analyzing angle-resolved data')
		initial_guesses, init_units = get_initial_from_yaml(config_params)
		angle_data, absorbance_data = get_angle_data_from_dir(spectral_data)
		sample, params = get_sample_params(spectral_data)
		write_angle_spec_to_file(angle_data, sample, output_path)
		angle, Elp, Eup, E_vib, new_wavenum, new_intens, resid = lorentzian_parsing(angle_data, absorbance_data, fit_func, bounds)
		#TODO: Also return the error
		splitting_fit = splitting_least_squares(initial_guesses, angle, Elp, Eup)

		E_0 = splitting_fit.x[0]
		E_vib = splitting_fit.x[1]
		Rabi = splitting_fit.x[2]
		refractive_index = splitting_fit.x[3]
		
# 		print(splitting_fit)
		print("Initial guesses (converted from {}): \n".format(init_units), initial_guesses)
		print("Fitting results:")
		print("E_0 =", E_0)
		print("E_vib =", E_vib)
		print("Rabi =", Rabi)
		print("n =", refractive_index)
		
		rad = [a * np.pi/180 for a in angle]  # convert degrees to radians
		E_cav = cavity_mode_energy(rad, E_0, refractive_index)  # calculate cavity mode
		
		write_dispersion_to_file(angle, Eup, Elp, E_vib, E_cav, sample, output_path)
		write_splitting_fit_to_file(splitting_fit.x, sample, init_units, output_path)
		
		for i, d in enumerate(angle_data):
			
			angle = d[0]
			wavenum = new_wavenum[i]
			intens = new_intens[2]
			residual_data = resid[i]
			write_peakfit_residuals_to_file(wavenum, residual_data, angle, output_path)
			i+=1

	else:
		print('No input data found')
		sys.exit()


if __name__ == '__main__':
	main()
