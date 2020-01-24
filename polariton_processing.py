#! /anaconda3/bin/python

"""
Polariton Data Processing
Author: Garrek Stemo
Date Created: November 18, 2019
Date Updated: December 17, 2019
Description: This program takes spectral data in .csv format
and performs curve fitting and other analysis, including graphical
output.
"""

import argparse
import os
import re
import sys
import csv
import numpy as np
from scipy import optimize
from scipy.interpolate import interp1d
from scipy import constants
import matplotlib.pyplot as plt
from ruamel_yaml import YAML
import pdb

yaml = YAML()


class Lorentzian:
	def __init__(self, amplitude=None, x0=None, y0=None, gamma=None):
		"""These parameters function as initial guesses for the fit function."""
		if amplitude == None:
			self.amplitude = 0.5
		if x0 == None:
			self.x0 = 1000.
		if y0 == None:
			self.y0 = 0.
		if gamma == None:
			self.gamma = 30.

	def set_amplitude(self, new_amp):
		self.amplitude = new_amp

	def set_x0(self, new_x0):
		self.x0 = new_x0

	def set_y0(self, new_y0):
		self.y0 = new_y0

	def set_gamma(self, new_gamma):
		self.gamma = new_gamma

	def lor_func(self, x):
		"""Lorentzian function for numpy array, x"""
		lor = self.y0 + self.amplitude * self.gamma**2 / ((x - self.x0)**2 + self.gamma**2)
		return lor


# ========== Unit conversions ========== #

def deg_to_rad(angles):
	"""Convert degrees to radians."""
	angles = [a * np.pi/180 for a in angles]
	return angles

def wavenum_to_wavelen(wavenum):
	"""cm^-1 to micrometers"""
	wavelen = (1/wavenum) * 10000
	return wavelen
	
def joule_to_ev(joule):
	ev = joule / constants.elementary_charge
	return ev
	
def wavenum_to_joule(wavenum):
	"""cm^-1 to photon energy"""
	cm_to_m = 1/100
	joules = constants.h * constants.c * (wavenum / cm_to_m)
	return joules

def wavenum_to_ev(wavenum):
	"""cm^-1 to eV units"""
	energy = wavenum_to_joule(wavenum)
	ev = joule_to_ev(energy)
	return ev

# ========== Get and write params and data ========== #

def get_bounds_from_yaml(yaml_config):
	"""Takes a yaml config file and gets the upper and lower bound
		used for truncating and Lorentzian fitting later."""
	with open(yaml_config, 'r') as yml:
		config = yaml.load(yml)
		
	lower = config['bounds']['lower']
	upper = config['bounds']['upper']
	return [lower, upper]
	
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
	
	return [E_0_guess, E_vib_guess, Rabi_guess, n_guess], initial_units
	

def get_data(spectral_data):
	"""Retrieve csv spectral data and store in array."""

	k = []
	I = []
	with open(spectral_data, 'r', errors='ignore') as spectrum:
		csvreader = csv.reader(spectrum)
		num_header_rows = 19   # Yeah this is hard coded. Not so great.
		for s in range(num_header_rows):
			next(csvreader, None)
		for row in csvreader:
			if row:
				wavenum = float(row[0])
				Inten = float(row[1])
				k.append(wavenum)
				I.append(Inten)
			else:
				break

	return np.array(k), np.array(I)

def get_angle_data_from_dir(directory):
	"""Extracts angle-resolved and absorbance data from each file in the supplied directory"""

	abs_k = [] 			 # absorbance wavenumbers
	abs_I = []  		 # absorbance intensities
	angle_data = []  # List with structure [[angle, [wavenumbers, intensities]], [...], ...]

	deg_str = 'deg'
	abs_str = 'Abs'

	for spectrum in os.listdir(directory):
		if spectrum.endswith('.csv'):
			spec_file = str(directory) + '/' + str(spectrum)
			if deg_str in spectrum:
				# Get the angle of the measurement from file string
				deg_s = spectrum.find(deg_str) + len(deg_str)
				deg_e = spectrum.find('_', deg_s)
				deg = int(spectrum[deg_s:deg_e])
				wavenum, intensity = get_data(spec_file)
				angle_data.append([deg, wavenum, intensity])

			# Get absorbance data if the file exists (it should always exist)
			if abs_str in spectrum:
				abs_k, abs_I = get_data(spec_file)

	absorbance_data = [abs_k, abs_I]
	angle_data.sort()
<<<<<<< HEAD
	
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
=======
>>>>>>> cavity_sim_validation

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

def truncate_data(xdata, ydata, bound1, bound2):
	"""Truncate data to isolate desired peaks for fitting"""

	i = 0
	MAX = len(xdata)
	lower_pt = 0
	upper_pt = 0

	while i < MAX:
		if np.round(xdata[i], 0) == float(bound1):
			lower_pt = i
		if np.round(xdata[i], 0) == float(bound2):
			upper_pt = i
		i+=1

	return xdata[lower_pt:upper_pt], ydata[lower_pt:upper_pt]

def write_angle_spec_to_file(angle_data_list, sample_name):
	"""Takes angles list, wavenumber list, and intensity list.
		Writes angle-resolved data to csv file with
		wavenumber in leftmost column (x-axis), and intensities
		for each degree (y-axes)."""

	angles = []
	wavenumbers = angle_data_list[0][1]
	num_wavenums = len(wavenumbers)

	angle_res_file = sample_name + '_angle-resolved_spectra.csv'
	output = os.path.join(os.path.abspath(args.output), angle_res_file)
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
<<<<<<< HEAD
=======

>>>>>>> cavity_sim_validation
			for item in angle_data_list:
				# Cycle through data point for each degree, make a row
				row.append(item[2][i])

			filewriter.writerow(row)
			i+=1

	print('')
	print('Wrote angle-resolved spectra results to {}\n'.format(output))

	return 0

def write_dispersion_to_file(angles, E_up, E_lp, E_vib, E_cav, sample_name):
	"""Takes angles and energies (wavenumbers) for upper and lower
	   polaritons, vibrational energy, calculated cavity mode
	   and writes all that data to a csv file.
	   REMEMBER: change plot_dispersion in plots.py if output file format changes."""

	dispersion_file = sample_name + '_dispersion.csv'
	output = os.path.join(os.path.abspath(args.output), dispersion_file)
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
	print('')
	print('Wrote dispersion results to {}'.format(output))
	return 0

def write_splitting_fit_to_file(least_squares_results, sample_name, units):
	"""Takes results of nonlinear least squares fit for Rabi splitting
	   (E_cav_0, Rabi, n, E_vib) and writes them to a file for 
	   later plotting and whatnot."""
	
	E_cav_0, E_vib, Rabi, n_ = least_squares_results
	least_squares_file = sample_name + '_splitting_fit.csv'
	output = os.path.join(os.path.abspath(args.output), least_squares_file)
	
	with open(output, 'w') as f:
		filewriter = csv.writer(f, delimiter=',')
		filewriter.writerow(['units', units])  #TODO: What units to write?
		filewriter.writerow(['E_cav_0', E_cav_0])
		filewriter.writerow(['Rabi', Rabi])
		filewriter.writerow(['n', n_])
		filewriter.writerow(['E_vib', E_vib])

	print('Wrote results of nonlinear least squares to {}'.format(output))
	return 0


# ========== Fitting Functions ========== #


def lor_1peak(x, A, x0, gamma, y0):
	"""Lorentzian fitting function for a single peak.
	   This is needed separate from the function in the Lorentzian class
	   in order to perform fits."""

	lor = y0 + A * gamma**2 / ((x - x0)**2 + gamma**2)
	return lor

def lor_2peak(x, A1, x0_1, gamma1, y0_1, A2, x0_2, gamma2, y0_2):
	"""Lorentzian fitting function for two peaks"""
	lor2 = y0_1 + A1 * gamma1**2 / ((x - x0_1)**2 + gamma1**2) \
		   + y0_2 + A2 * gamma2**2 / ((x - x0_2)**2 + gamma2**2)
	return lor2

def coupled_energies(theta, E0, Ee, Rabi, n_eff, branch=0):
	
	Ec = E0 / np.sqrt(1 - (np.sin(theta) /n_eff)**2)
	
	if branch == 0:
		E_coupled = 0.5*(Ee + Ec) - 0.5*np.sqrt(4*Rabi**2 + (Ee - Ec)**2)
		
	elif branch == 1:
		E_coupled = 0.5*(Ee + Ec) + 0.5*np.sqrt(4*Rabi**2 + (Ee - Ec)**2)
	
	return E_coupled

def cavity_mode_energy(angle, E_0, n_eff):
	"""E0 is some initial (0 deg incidence) cavity mode energy, angle is in radians,
	   n_eff is the effective refractive index of the material."""
	E_cav = E_0 / np.sqrt(1 - (np.sin(angle) /n_eff)**2)
	return E_cav

def error_f(x, theta, Elp_data, Eup_data):
	"""Here, we write the eigenvalues of the energy Hamiltonian for 
	   coupled cavity-vibration system. We minimize this with experimental data."""

	En = coupled_energies(theta, *x, branch=0)
	Ep = coupled_energies(theta, *x, branch=1)
	
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
	
	E_0, Rabi, n_eff, E_vib = x
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


# ========== Fitting Procedures ========== #

def absorbance_fitting(wavenum, intensity, bounds):
	"""Takes list of wavenumber and intensity lists for x, y axes.
	   Takes bounds to truncate for fitting.
	   Returns Lorentzian function class with fitted parameters for absorbance."""

	low_bound = bounds[0]
	up_bound = bounds[1]
	k, I = truncate_data(wavenum, intensity, low_bound, up_bound)
	lor = Lorentzian()
	center = low_bound + 1/2 * (up_bound - low_bound)
	lor.set_x0(center)
	
	p0=[lor.amplitude, lor.x0, lor.gamma, lor.y0]

	popt, pconv = optimize.curve_fit(lor_1peak, k, I, p0)

	lor.amplitude = popt[0]
	lor.x0 = popt[1]
	lor.gamma = popt[2]
	lor.y0 = popt[3]
	lor.lor_func(k)

	return lor

def polariton_fitting(wavenum, intensity, lorz1, lorz2):
	"""Fit the curve with data and Lorentzian class"""

	fitting_func = lorz1.lor_func(wavenum)		

# 	amp1 = lorz1.amplitude
# 	x01 = lorz1.x0
# 	g1 = lorz1.gamma
# 	y01 = lorz1.y0
	peak1_err = []

# 	amp2 = lorz2.amplitude
# 	x02 = lorz2.x0
# 	g2 = lorz2.gamma
# 	y02 = lorz2.y0
	peak2_err = []

	p0 = [lorz1.amplitude, lorz1.x0, lorz1.gamma, lorz1.y0,
		  lorz2.amplitude, lorz2.x0, lorz2.gamma, lorz2.y0]

	lor_fit = []

	try:
		popt, pcov = optimize.curve_fit(lor_2peak, wavenum, intensity, p0)

		lorz1.amplitude, lorz1.x0, lorz1.gamma, lorz1.y0 = popt[0:4]
		peak1_args = popt[0:4]
		peak2_args = popt[4:8]
		lorz2.amplitude, lorz2.x0, lorz2.gamma, lorz2.y0  = popt[4:8]
# 		print(np.diag(pcov))
		err = np.sqrt(np.diag(pcov))
		peak1_err = err[0:4]
		peak2_err = err[4:8]

		lor_fit = lor_2peak(wavenum, *peak1_args, *peak2_args)

	except RuntimeError:
		print("Something went wrong with finding a fit (RuntimeError).")
		return 1

	amp1 = np.round([lorz1.amplitude, peak1_err[0]], 2)
	center1 = np.round([lorz1.x0, peak1_err[1]], 2)
	gamma1 = np.round([lorz1.gamma, peak1_err[2]], 2)

	amp2 = np.round([lorz2.amplitude, peak2_err[0]], 2)
	center2 = np.round([lorz2.x0, peak2_err[1]], 2)
	gamma2 = np.round([lorz2.gamma, peak2_err[2]], 2)

# 	print('')
# 	print('-'*10, 'Lower Polariton', '-'*10)
# 	print("Amplitude = ", amp1[0], u'\u00b1', amp1[1])
# 	print("Center = ", center1[0], u'\u00b1', center1[1])
# 	print('sigma = ', gamma1[0], u'\u00b1', gamma1[1])
#
# 	print('')
# 	print('-'*10, 'Upper Polariton', '-'*10)
# 	print("Amplitude = ", amp2[0], u'\u00b1', amp2[1])
# 	print("Center = ", center2[0], u'\u00b1', center2[1])
# 	print('sigma = ', gamma2[0], u'\u00b1', gamma2[1])

	return lor_fit, lorz1, lorz2


def lorentzian_parsing(angle_data, absor_data, bounds):
	"""Takes raw angle-resolved spectra, absorbance spectrum, list with upper, lower bound.
	   Returns lists of angles, upper/lower polariton peak positions (cm-1),
	   absorbance peak position. This info is used for dispersion curves later."""

	low_bound = bounds[0]
	up_bound = bounds[1]

	# Initialize classes to store Lorentzian fit data
	lorz1 = Lorentzian()
	lorz2 = Lorentzian()

	# Initial polariton position guesses
	x01 = low_bound + 1/3 * (up_bound - low_bound)
	x02 = low_bound + 2/3 * (up_bound - low_bound)
	lorz1.set_x0(x01)
	lorz2.set_x0(x02)

	angles = [] 	  				# List of angles from experiment
	wavenumbers = angle_data[0][1]  # List of list of wavenumbers from experiment
	intensities = [] 				# List of list of intensities from experiment
	lor_fits = []  	  				# List of Lorentzian fits for each data set
	lower_pol = []	  				# List of Lorentzian classes for lower polariton
	upper_pol = []	  				# List of Lorentzian classes for upper polariton
	
	abs_k = absor_data[0]
	abs_I = absor_data[1]
	abs_lor = absorbance_fitting(abs_k, abs_I, bounds)

	for d in angle_data:

		angles.append(d[0])
		wavenum = d[1]
		intensity = d[2]

		wavenum, intensity = truncate_data(wavenum, intensity, low_bound, up_bound)

		# Fit the data to find upper and lower polaritons.
		fit, lor_lower, lor_upper = polariton_fitting(wavenum, intensity, lorz1, lorz2)
		lorz1 = lor_lower
		lorz2 = lor_upper

# 		wavenumbers.append(wavenum)
		intensities.append(intensity)
		lor_fits.append(fit)

		# For now, just using the x0 position of the peak
		lower_pol.append(lor_lower.x0)
		upper_pol.append(lor_upper.x0)

	return angles, lower_pol, upper_pol, abs_lor.x0


def splitting_least_squares(initial, angles, Elp, Eup):
	"""Takes initial guesses:[E_cav_0, E_vib, Rabi, n].
	   Takes angles (in degrees), and experimental upper and lower polariton data.
	   Returns nonlinear least squares fit."""

	angles = deg_to_rad(angles)
	#TODO: What units do we really want here? Unitless to convert later?
	Elp = [wavenum_to_ev(i) for i in Elp]
	Eup = [wavenum_to_ev(i) for i in Eup]

# 	vib_lb = Elp[-1]
# 	vib_ub = Eup[0]
# 	E0_lb = Elp[0]
# 	E0_ub = Eup[0]
# 	lb = [E0_lb, 0., 1.5, vib_lb]
# 	ub = [E0_ub, 0.5, 2.0, vib_ub]
# 	bound_vals = (lb, ub)
	print("Performing default nonlinear least squares fit.")
	optim = optimize.least_squares(error_f,
								   x0=initial,
								   args=(angles, Elp, Eup))
	return optim


def cavity_modes(bounds):
	"""Use at your own peril."""
	degree = []  # degree paired with spectrum file
	cavity_dir = args.cavity_mode
	for spectrum in os.listdir(cavity_dir):
		if spectrum.endswith('.csv'):
			spec_file = str(cavity_dir) + '/' + str(spectrum)
			deg_str = 'deg'
			if deg_str in spectrum:
				# Get the angle of the measurement from file string
				deg_s = spectrum.find(deg_str) + len(deg_str)
				deg_e = spectrum.find('_', deg_s)
				deg = float(spectrum[deg_s:deg_e])

				degree.append([deg, spec_file])
	degree.sort()

	low_bound = bounds[0]
	up_bound = bounds[1]
	x0 = float(args.cav_center)

	lor = Lorentzian()
	lor.set_x0(x0)
	angles = [d[0] for d in degree]
	wavenumbers = []  # List of list of wavenumbers
	intensities = []  # List of list of intensities
	lor_fits = []  	  # List of Lorentzian fits for each data set
	modes = []		  # List of Lorentzian classes for cavity modes

	for d in degree:
		spectrum_file = d[1]

		wavenum, intensity = get_data(spectrum_file)

		# --- These are for simulated data from Pistachio --- #
# 		wavenum = [np.power(10, 4)/w for w in wavenum]
# 		wavenum = list(reversed(wavenum))
# 		intensity = list(reversed(intensity))
		# --------------------------------------------------- #

		wavenum, intensity = truncate_data(wavenum, intensity, low_bound, up_bound)

		try:
			popt, pconv = optimize.curve_fit(lor_1peak, wavenum, intensity,
											 p0=[lor.amplitude, lor.x0, lor.gamma, lor.y0])
			lor.amplitude, lor.x0, lor.gamma, lor.y0 = popt

			lor_fit = lor_1peak(wavenum, lor.amplitude, lor.x0, lor.gamma, lor.y0)
			wavenumbers.append(wavenum)
			intensities.append(intensity)
			lor_fits.append(lor_fit)
			modes.append(lor.x0)

		except RuntimeError:
			print("Cavity mode runtime error.")
			return 1

	return angles, modes


# ====== Plotting (might move to separate module) ====== #

def plot_absorbance(k, I, lor, bounds):

	low_bound = bounds[0]
	up_bound = bounds[1]

	fig, (ax, axf) = plt.subplots(2)

	ax.plot(k, I)
	ax.set_xlim([int(up_bound), int(low_bound)])

	axf.plot(k, lor.lor_func(k))
	axf.set_xlim([int(up_bound), int(low_bound)])

	ax.set_xlabel(r'Wavenumber (cm$^{-1}$)')
	ax.set_ylabel('Transmission (%)')
	axf.set_xlabel(r'Wavenumber (cm$^{-1}$)')
	axf.set_ylabel('Transmission (%)')

	plt.show()

	return 0

def plot_polariton(x_, y_, fit_func):
	"""Makes a single plot for polariton data"""

	fig, ax = plt.subplots()
	ax.plot(x_, y_)
	ax.plot(x_, fit_func)

def main():

	spectral_data = args.spectral_data
	config_params = args.config
	bounds = get_bounds_from_yaml(config_params)
	bounds.sort()
	


	if args.polariton:
		print('Fitting double-peak Lorentzian')
		x, y, fit = polariton_fitting(spectral_data)
		plot_polariton(x, y, fit)

	elif args.absorbance:
		print('Fitting single-peak Lorentzian')
		k, I = get_data(spectral_data)
		lor = absorbance_fitting(k, I)
		plot_absorbance(k, I, lor)

	elif args.angleres:
		print('Analyzing angle-resolved data')
		initial_guesses, init_units = get_initial_from_yaml(config_params)
		ang_data, abs_data = get_angle_data_from_dir(spectral_data)
		sample, params = get_sample_params(spectral_data)

# 		write_angle_spec_to_file(ang_data, sample)
		ang, Elp, Eup, E_vib = lorentzian_parsing(ang_data, abs_data, bounds)
		#TODO: Also return the error
		splitting_fit = splitting_least_squares(initial_guesses, ang, Elp, Eup)

		E_0 = splitting_fit.x[0]
		E_vib = splitting_fit.x[1]
		Rabi = splitting_fit.x[2]
		refractive_index = splitting_fit.x[3]
		
# 		print(splitting_fit)
		print("Initial guesses ({}): ".format(init_units), initial_guesses)
		print("Fitting results:")
		print("E_0 =", E_0)
		print("E_vib =", E_vib)
		print("Rabi =", Rabi)
		print("n =", refractive_index)

		
		rad = [a * np.pi/180 for a in ang]  # convert degrees to radians
		E_cav = cavity_mode_energy(rad, E_0, refractive_index)  # calculate cavity mode
		
		write_dispersion_to_file(ang, Eup, Elp, E_vib, E_cav, sample)
		write_splitting_fit_to_file(splitting_fit.x, sample, init_units)
		
	else:
		print('No input data found')
		sys.exit()


if __name__ == '__main__':

	parser = argparse.ArgumentParser()

	spectrum_help = "csv file or directory containing spectral information."
	output_help = "Path for output data directory."
	config_help = "Yaml file to set Lorentz fit bounds and \
					least squares fit initial guesses for \
					0 degree incidence cavity mode energy, \
					Rabi splitting value, \
					vibrational excitation energy."
	cavity_help = "csv file containing angle-tuned cavity mode data."
	cav_cen_help = "Initial guess for x position of mode center."
	spec_dir_help = "Directory of angle-resolved polariton spectral data."
	pol_help = "Boolean. Lorentzian fit for a single spectrum file."
	abs_help = "Boolean indicating absorbance single peak data."
	angle_help = "Boolean indicating directory contains angle-resolved data."
	

	parser.add_argument('spectral_data', help=spectrum_help)
	parser.add_argument('output', help=output_help)
	parser.add_argument('-CF', '--config', help=config_help)
	parser.add_argument('-C', '--cavity_mode', help=cavity_help)
	parser.add_argument('-E', '--cav_center', help=cav_cen_help)
	parser.add_argument('-P', '--polariton', action='store_true', help=pol_help)
	parser.add_argument('-A', '--absorbance', action='store_true', help=abs_help)
	parser.add_argument('-T', '--angleres', action='store_true', help=angle_help)

	args = parser.parse_args()

	main()
