"""
Polariton Data Processing
Author: Garrek Stemo
Date Created: November 18, 2019
Date Updated: November 27, 2019
Description: This program takes spectral data in .csv format
and performs curve fitting and other analysis, including graphical
output.
"""

import argparse
import os
import sys
import csv
import numpy as np
from scipy import optimize
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


class Lorentzian:
	def __init__(self, amplitude=None, x0=None, y0=None, gamma=None):
		"""These parameters function as initial guesses for the fit function."""
		if amplitude == None:
			self.amplitude = 1.
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
	

def get_data(spectral_data):
	"""Retrieve csv spectral data and store in array."""
	
	k = []
	I = []
	with open(spectral_data, 'r', errors='ignore') as spectrum:
		csvreader = csv.reader(spectrum)
		for s in range(19):
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
	

def polariton_fitting(wavenum, intensity, lorz1, lorz2):
	"""Fit the curve with data and Lorentzian class"""

	fitting_func = lorz1.lor_func(wavenum)
	
	amp1 = lorz1.amplitude
	x01 = lorz1.x0
	g1 = lorz1.gamma
	y01 = lorz1.y0
	peak1_err = []
	
	amp2 = lorz2.amplitude
	x02 = lorz2.x0
	g2 = lorz2.gamma
	y02 = lorz2.y0
	peak2_err = []

	lor_fit = []

	try:
		popt, pconv = optimize.curve_fit(lor_2peak, wavenum, intensity,
										 p0=[amp1, x01, g1, y01, amp2, x02, g2, y02])
		
		lorz1.amplitude, lorz1.x0, lorz1.gamma, lorz1.y0 = popt[0:4]
		peak1_args = popt[0:4]
		peak2_args = popt[4:8]
		lorz2.amplitude, lorz2.x0, lorz2.gamma, lorz2.y0  = popt[4:8]
		err = np.sqrt(np.diag(pconv))
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


def plot_polariton(x_, y_, fit_func):
	"""Makes a single plot for polariton data"""
	
	fig, ax = plt.subplots()
	ax.plot(x_, y_)
	ax.plot(x_, fit_func)
	

def angle_data():
	"""Input: directory of angle-resolved spectra with polariton double peaks"""

	degree = []  # degree paired with spectrum file
	directory = args.spectral_data
	abs_k = []
	abs_I = []
	abs_lor = None # wavenum, intensity, Lorentzian class
	for spectrum in os.listdir(directory):
		if spectrum.endswith('.csv'):
			spec_file = str(directory) + '/' + str(spectrum)
			if 'deg' in spectrum:
				# Get the angle of the measurement from file string
				deg_s = spectrum.find('deg') + 3
				deg_e = spectrum.find('_', deg_s)
				deg = int(spectrum[deg_s:deg_e])
				degree.append([deg, spec_file])
				
			if 'Abs' in spectrum:
				# Get absorbance data if the file exists
				 abs_k, abs_I, abs_lor = absorbance_fitting(spec_file)
	degree.sort()

	up_bound = float(args.upperbound)
	low_bound = float(args.lowerbound)
	x01 = low_bound + 1/3 * (up_bound - low_bound)
	x02 = low_bound + 2/3 * (up_bound - low_bound)
	
	# Make classes to store Lorentzian fit data
	lorz1 = Lorentzian()
	lorz2 = Lorentzian()
	
	lorz1.set_x0(x01)
	lorz2.set_x0(x02)
	
	angles = [int(d[0]) for d in degree]
	wavenumbers = []  # List of list of wavenumbers
	intensities = []  # List of list of intensities
	lor_fits = []  	  # List of Lorentzian fits for each data set
	upper_pol = []	  # List of Lorentzian classes for upper polariton
	lower_pol = []	  # list of Lorentzian classes for lower polariton
		
	for d in degree:
		spectrum_file = d[1]
# 		print('='*40)
# 		print(d[0], 'degrees')
# 		print(spectrum_file)
		
		wavenum, intensity = get_data(spectrum_file)
		wavenum, intensity = truncate_data(wavenum, intensity, low_bound, up_bound)

		fit, lor_upper, lor_lower = polariton_fitting(wavenum, intensity, lorz1, lorz2)
		lorz1 = lor_lower
		lorz2 = lor_upper
		wavenumbers.append(wavenum)
		intensities.append(intensity)
		lor_fits.append(fit)
		
		# For now, just using the x0 position of the peak
		upper_pol.append(lor_upper.x0)
		lower_pol.append(lor_lower.x0)
	
# 	fig0, (ax, axf) = plt.subplots(2)
	fig0, ax = plt.subplots()
	
	for i in range(len(wavenumbers)):
		ax.plot(wavenumbers[i], intensities[i])
		ax.set_xlim([int(args.upperbound), int(args.lowerbound)])
		ax.set_ylim(0., 0.9)

# 		axf.plot(wavenumbers[i], lor_fits[i])
# 		axf.set_xlim([int(args.upperbound), int(args.lowerbound)])
# 		axf.set_ylim(0., 0.9)

	
	ax.set_ylabel('Transmission (%)', fontsize=12)
	ax.set_xlabel(r'Wavenumber (cm$^{-1}$)', fontsize=12)
# 	axf.set_ylabel('Transmission (%)', fontsize=12)
	ax.set_title('Angle-tuned measurements for Neat DPPA', fontsize=14)
	
	fig1, ax = plt.subplots()
# 	ax.scatter(angles, upper_pol, s=10, c='black')
# 	ax.scatter(angles, lower_pol, s=10, c='black')
	ax.set_xticks(angles)
# 	ax.set_ylim(0., 0.8)
	ax.set_ylim(2060, 2160)
	
	abs_pos = abs_lor.x0
	abs_y = np.zeros(len(angles))
	abs_y = [abs_pos for i in abs_y]
# 	ax.plot(angles, abs_y, 
# 			linestyle='--', c='black', label='N=N=N stretch mode')
	
# 	cav_angles, cav_modes = cavity_modes()
# 	ax.plot(cav_angles, cav_modes, 
# 			linestyle='--', c='blue')
	
	ax.set_title('Angle-tuned empty cavity', fontsize=14)
	ax.set_xlabel('Angle (degrees)', fontsize=12)
	ax.set_ylabel(r'Wavenumber (cm$^{-1}$)', fontsize=12)
# 	ax.legend()

	# Save figures as PDFs	
# 	fig0.savefig(os.path.join('/Users/garrek/Documents/project_graphics/Nov 28 Big Group Meeting/',
# 							  'angle-tuned_Neat_DPPA.pdf'),
# 							  bbox_inches='tight', 
# 							  dpi=300)
# 	fig1.savefig(os.path.join('/Users/garrek/Documents/project_graphics/Nov 28 Big Group Meeting/', 
# 							  'dispersion_vacant_cavity.pdf'), dpi=300)
# 	plt.show()


def Ecav(A, B):
	
	k = np.linspace(0, np.pi/4, np.pi/100)
	
	return A + B*k**2

def cavity_modes_sim():
	"""Not Implemented"""
	
	Ep = 0.5*(Ecav + Evib) + np.sqrt(V**2 + (Evib - Ecav)**2)
	Em = 0.5*(Ecav + Evib) - np.sqrt(V**2 + (Evib - Ecav)**2)
	
	return 0

def interpolate_polariton():
	"""Not Implemented"""


def cavity_modes():
	"""Use at your own peril."""
	degree = []  # degree paired with spectrum file
	cavity_dir = args.cavity_mode
	for spectrum in os.listdir(cavity_dir):
		if spectrum.endswith('.csv'):
			spec_file = str(cavity_dir) + '/' + str(spectrum)
			if 'deg' in spectrum:
				# Get the angle of the measurement from file string
				deg_s = spectrum.find('deg') + 3
				deg_e = spectrum.find('_', deg_s)
				deg = float(spectrum[deg_s:deg_e])

				degree.append([deg, spec_file])	
	degree.sort()

	up_bound = float(args.upperbound)
	low_bound = float(args.lowerbound)
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
# 		print('='*40)
# 		print(d[0], 'degrees')
# 		print(spectrum_file)
		
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

def absorbance_fitting(spectrum_file):
	"""Fits and plots absorbance data
	   Returns wavenumbers, intensities, and Lorentzian function class with
	   fitted parameters."""
	
	k, I = get_data(spectrum_file)
	up_bound = float(args.upperbound)
	low_bound = float(args.lowerbound)
	k, I = truncate_data(k, I, low_bound, up_bound)
	lor = Lorentzian()
	center = low_bound + 1/2 * (up_bound - low_bound)
	lor.set_x0(center)
	
	popt, pconv = optimize.curve_fit(lor_1peak, k, I,
								 p0=[lor.amplitude, lor.x0, lor.gamma, lor.y0])
	
	lor.amplitude = popt[0]
	lor.x0 = popt[1]
	lor.gamma = popt[2]
	lor.y0 = popt[3]
	lor.lor_func(k)
	print("Vibrational mode: ", lor.x0)

	return k, I, lor
	
	
def plot_absorbance(k, I, lor):

	up_bound = float(args.upperbound)
	low_bound = float(args.lowerbound)

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


def get_params(file_name):
	"""Gets spectrum parameters from file name"""


def main():
	
	if args.polariton:
		print('Fitting double-peak Lorentzian')
		x, y, fit = polariton_fitting(args.spectral_data)
		plot_polariton(x, y, fit)
		
	elif args.absorbance:
		print('Fitting single-peak Lorentzian')
		k, I, lor = absorbance_fitting(args.spectral_data)
		plot_absorbance(k, I, lor)
		
	elif args.angleres:
		print('Fitting angle-resolved data')
		angle_data()
		
	else:
		print('No input data found')
		sys.exit()


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	
	spectrum_help = "csv file containing spectral information."
	cavity_help = "csv file containing angle-tuned cavity mode data."
	cav_cen_help = "Initial guess for x position of mode center."
	spec_dir_help = "Directory of angle-resolved polariton spectral data."
	pol_help = "Boolean. Lorentzian fit for a single spectrum file."
	abs_help = "Boolean indicating absorbance single peak data."
	angle_help = "Boolean indicating directory contains angle-resolved data."
	lbound_help = "Lower x-axis bound for fitting."
	ubound_help = "Upper x-axis bound for fitting."
	
	parser.add_argument('spectral_data', help=spectrum_help)
	parser.add_argument('-C', '--cavity_mode', help=cavity_help)
	parser.add_argument('-CC', '--cav_center', help=cav_cen_help)
	parser.add_argument('-P', '--polariton', action='store_true', help=pol_help)
	parser.add_argument('-A', '--absorbance', action='store_true', help=abs_help)
	parser.add_argument('-T', '--angleres', action='store_true', help=angle_help)
	parser.add_argument('-L', '--lowerbound', help=lbound_help)
	parser.add_argument('-U', '--upperbound', help=ubound_help)
	
	args = parser.parse_args()
	
	main()
