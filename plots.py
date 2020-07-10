#! /anaconda3/bin/python

import polariton_processing as pp
import argparse
import csv
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from scipy import constants
from scipy import stats
import convert_unit


# ================ Get data from files ================ #

def get_splitting_results(splitting_file):
	"""These parameters were generated by nonlinear least squares fitting."""
	params = {}
	with open(splitting_file, 'r') as sf:
		csvreader = csv.reader(sf)
		for row in csvreader:
			try:
				params[row[0]] = float(row[1])
			except ValueError:
				params[row[0]] = row[1]

	return params
	
def get_concentration_data(concentration_file):
	"""Get sorted data generated from concentration_analysis"""
	conc = []
	rabi = []
	with open(concentration_file, encoding='utf-8', newline='',) as f:
		csvreader = csv.reader(f)
		next(csvreader, None)
		for row in csvreader:
			conc.append(float(row[0]))
			rabi.append(float(row[1]))
	
	return conc, rabi

# =================================================================== #
# =================      Plotting functions 	 ==================== #
# =================================================================== #

def tmm_plots(sim_path, save_plot=None):
	"""
	Plot transmission from transfer matrix simulation.
	"""
	#TODO: User specifies transmittance, reflectance, absorbance.
	#TODO: User specifies center wavenumber

	fig, ax = plt.subplots(figsize=(20, 10))
	gs1 = gridspec.GridSpec(3, 1)
	gs1.update(wspace=0.025, hspace=0.005)
		
	sim_file_list = os.listdir(sim_path)
	for sim in sim_file_list:
	
		wavelength = []
		wavenumber = []
		transmittance = []
		reflectance = []
		absorptance = []
		field = []
		
		sim_file = os.path.join(sim_path, sim)
		with open(sim_file, 'r', encoding='utf8', newline='') as f:
			reader = csv.reader(f)
			next(reader, None)
			for row in reader:
				wavelength.append(float(row[0]))
				transmittance.append(float(row[1]))
				reflectance.append(float(row[2]))
				absorptance.append(float(row[3]))

		wavenumber = [10**4/wl for wl in wavelength]
		ax.plot(wavenumber, transmittance,
				color='k',
				linewidth=1,
				label="transfer matrix")
# 	if args.transmission:
# 		ax.plot(wl_T_data, T_data, linestyle="dashed", color='#FF5733', label="downloaded data")
	ax.set_ylabel('Transmittance %', fontsize=16)
	ax.set_xlabel(r'Wavenumber (cm$^{-1}$)', fontsize=16)
	ax.tick_params(axis='both', which='both', direction='in', right=True, top=True, labelsize=14)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
	ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
	ax.set_xlim(1500, 500)
	ax.set_ylim(0.0, 5.1)
	plot_title = "Transfer Matrix Method"
	plt.suptitle(plot_title, fontsize=24)
# 	plt.subplots_adjust(top=1.0)

	if save_plot:
		print("Saving figure as pdf to {}".format(save_plot))
		fig.savefig(save_plot, bbox_inches='tight')

	else:
		plt.show()


def reference_data(data_file):
	"""Gets reference data downloaded from
	websites. Filmetrics.com data are in nanometers"""
	wavelength = []
	Y = []
	unit = 1  # sets order of magnitude (nm or um)
	with open(data_file, 'r') as ref:
	
		reader = None
		if 'filmetrics' in str(ref):
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


# ====== Plot Polariton Data and Dispersion Curves ====== #

def plot_spectra(file_prefix, spectra_file, excitation=None, save_dir=None):
	"""Takes csv file with spectral data produced by
	   write_angle_spec_to_file or write_dispersion_to_file functions in
	   polariton_processing module"""
	
	wavenumbers = []
	intensities = []
	angles = []

	with open(spectra_file) as sfile:
		csvreader = csv.reader(sfile)
		header = next(csvreader)
		deg_str = 'deg'
		
		for deg in header[1:]:
			ang_start = deg.find(deg_str) + len(deg_str)
			ang_end = deg.find(' ', ang_start)
			ang = int(deg[ang_start:ang_end])
			angles.append(ang)
			intensities.append([ang, []])

		for row in csvreader:
			wavenumbers.append(float(row[0]))
			
			inten_data = row[1:]
			for idx, a in enumerate(intensities):
				a[1].append(float(inten_data[idx]))

	fig, ax = plt.subplots()
	
	y_offset = 0.
	i = 0
	MAX = len(angles) - 1
	theta1 = str(angles[0])
	theta2 = str(angles[MAX])
	while i < len(angles):

		y_vals = intensities[i][1]
		y_vals = [y+y_offset for y in y_vals]

		deg_label = str(angles[i])
		ax.plot(wavenumbers, y_vals,
				color='black',
				linewidth=0.5,
				linestyle='dashed',
				label=deg_label)
		y_offset += 0.
# 		i-=1
		i+=1
		
	if excitation:
		ax.axvline(x=excitation)

	# Figure formatting
	ax.tick_params(axis='both', which='both', direction='in', right=True, top=True)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
	ax.yaxis.set_minor_locator(ticker. AutoMinorLocator(5))
	ax.set_xlim([3000, 1400])
	ax.set_ylim([0, 2])

	# Annotate	
	#TODO: Don't hard code these labels or xlim, ylim.
# 	xy_pt1 = (2350, 0.09)
# 	xy_pt2 = (2350, 2.7)

	angle_info_str = "Angle range: ({}, {})".format(theta1, theta2)

# 	ax.annotate(r'$\theta$ = {}$\degree$'.format(theta1), xy=xy_pt1)
# 	ax.annotate(r'$\theta$ = {}$\degree$'.format(theta2), xy=xy_pt2)
	ax.text(0.45, 1.0, angle_info_str, 
			horizontalalignment='left',
			verticalalignment='top',
			transform=ax.transAxes,
			bbox=dict(boxstyle='square', facecolor='white'))

	
	ax.set_xlabel(r'Wavenumber (cm$^{-1}$)')
	ax.set_ylabel('Transmission %')

	# Save file as PDF or something
	if save_dir:
		file_name = file_prefix + '_cascade_plot.pdf'
		output_file = os.path.join(save_dir, file_name)
		fig.savefig(output_file, bbox_inches='tight')
		print("Saved cascade plot to {}".format(output_file))
	else:
		plt.show()
	
	return 0
	

def plot_dispersion(file_prefix, dispersion_file, plot_units, splitting=None, save_dir=None):
	"""Takes dispersion curve file and plots wavenumber vs angle
	   for UP, LP, vibration mode, cavity mode"""
	
	if plot_units == 'ev':
		ylabel = 'Energy (eV)'
		unit_str = 'eV'
	elif plot_units == 'wl':
		ylabel = r'Wavelength ($\mu$m)'
		unit_str = 'wavelength'
	elif plot_units == 'wn':
		ylabel = r'Wavenumber (cm$^{-1}$)'
		unit_str = 'wavenumber'
	
	angles = []
	up = []
	lp = []
	vibration = []
	cavity = []
	
	with open(dispersion_file, 'r') as dfile:
		csvreader = csv.reader(dfile)
		next(csvreader)
		for row in csvreader:
			angles.append(int(row[0]))
			up.append(float(row[1]))
			lp.append(float(row[2]))
			vibration.append(float(row[3]))
			cavity.append(float(row[4]))
	current_units = 'wn' # wavenumber
	up = convert_unit.set_units(up, current_units, plot_units)
	lp = convert_unit.set_units(lp, current_units, plot_units)
	cavity = convert_unit.set_units(cavity, current_units, plot_units)
	vibration = convert_unit.set_units(vibration, current_units, plot_units)

	fig, ax = plt.subplots()

	mark_size = 10
	color1 = 'dimgray'
	color2 = 'black'
	color3 = 'deeppink'
	
	vib_label = 'molecular stretch'
	vib_xy = (0.89, vibration[0])
		
	cav_label = 'cavity dispersion'
	cav_xy = (0.89, cavity[-1])

	# Generate theoretical data
	n_points = 100
	t_min = -35
	t_max = 35
	theta_plot = np.linspace(t_min, t_max, n_points)
	theta_rad = [a * np.pi/180 for a in theta_plot]
	
	# Nonlinear least squares results labeling
	if splitting:
		
		# Get fitting params from text file and convert units
		params = get_splitting_results(splitting)
		E_cav0 = params['E_cav_0']
		n_eff = params['n']
		Rabi = params['Rabi']
		E_vib = params['E_vib']
		E_cav0, E_vib, Rabi = convert_unit.set_units([E_cav0, E_vib, Rabi], current_units, plot_units)
		
		# Generate curve from fitting data
		Ec = pp.cavity_mode_energy(theta_rad, E_cav0, n_eff)
		E_lp = pp.coupled_energies(theta_rad, E_cav0, E_vib, Rabi, n_eff, 0)
		E_up = pp.coupled_energies(theta_rad, E_cav0, E_vib, Rabi, n_eff, 1)

		e_vib_plot = np.full((n_points, ), E_vib)
		ax.plot(theta_plot, Ec, color=color2)
		ax.plot(theta_plot, E_up, color=color2)
		ax.plot(theta_plot, E_lp, color=color2)
		ax.plot(theta_plot, e_vib_plot, linestyle='dashed', color=color2)
		
		textstr = '\n'.join((
				  'Least Squares Fit',
				   r'$\Omega_R = %.2f$' % (Rabi),
				   r'$E_{cav,0} = %.1f$' % (E_cav0),
				   r'$E_{vib} = %.1f$' % (E_vib),
				   r'$n = %.2f$' % (n_eff)))
		if not save_dir:
			ax.text(0.45, 1.0, textstr, fontsize=10,
					horizontalalignment='left', verticalalignment='top',
					transform=ax.transAxes,
					bbox=dict(boxstyle='square', facecolor='white'))
	else:
		vibration = np.full((n_points, ), vibration[0])
		ax.plot(theta_plot, vibration,
			linestyle='dashed',
			color=color1,
			label=vib_label)

	#Plot experimental data
	ax.scatter(angles, up, s=mark_size, color=color3)
	ax.scatter(angles, lp, s=mark_size, color=color3)

	# Figure formatting
	ax.tick_params(axis='both', which='both', direction='in', right=True, top=True)
	angle_ticks = [-30, -20, -10, 0, 10, 20, 30]
	ax.set_xticks(angle_ticks, minor=True)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
	ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
	ax.set_xlabel(r'Incident angle (deg)')
	ax.set_ylabel(ylabel)

	# Some settings for saving to PDF
	if save_dir:
		ax.text(1.01, 0.5, textstr, fontsize=10,
				horizontalalignment='left', verticalalignment='top',
				transform=ax.transAxes,
				bbox=dict(boxstyle='square', facecolor='white'))
		file_name = file_prefix + '_dispersion_curve.pdf'
		output_file = os.path.join(save_dir, file_name)
		fig.savefig(output_file, bbox_inches='tight')
		print("Saved dispersion plot to {}".format(output_file))
	else:
		plt.show()
	
	return 0

def plot_splitting_concentration(concentration_file, new_units=None):
	"""Plots solute concentration vs vacuum Rabi splitting"""
	#TODO: Plot more than one pair of data
	concentrations, splittings = get_concentration_data(concentration_file)
	m, y0, rval, pval, std = stats.linregress(concentrations, splittings)

	c_fit = np.linspace(0., 1., 100)
	s_fit = m*c_fit + y0

# 	splittings = convert_unit.set_units(splittings, 'wn', 'ev')
	fig, ax = plt.subplots()
	
	ax.scatter(concentrations, splittings)
	ax.plot(c_fit, s_fit)
	
	# Set axes and ticks
	ax.set_xlim([0., 1.])
	ax.set_ylim(ymin=0.)
	ax.tick_params(axis='both', which='both', direction='in')
	ax.xaxis.set_major_locator(ticker.AutoLocator())
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
	ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
	
	ax.yaxis.set_major_locator(ticker.AutoLocator())
	ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
	ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
	
	ax.set_xlabel(r'$\sqrt{c}$')
	ax.set_ylabel(r'Splitting parameter (cm$^{-1}$)')
	
	# Annotate
	fit_str = '\n'.join(('Linear fit', '%.2f + %.2f x' % (m, y0)))
	bbox_props = dict(boxstyle='square', fc='white')
	ax.text(0.7, 10, fit_str, bbox=bbox_props)
	
	plt.show()

def plot_curve_fit_residuals(residuals_file):
	"""Plots residuals for spectral curve fitting"""

	x = []	
	residuals = []
	
	with open(residuals_file, 'r') as rf:
		csvreader = csv.reader(rf)
		for row in csvreader:
			x.append(float(row[0]))
			residuals.append(float(row[1]))
			
	fig, ax = plt.subplots()
	
	ax.scatter(x, residuals, s=1)
	ax.plot(x, np.zeros(len(x)), linestyle='dashed')
	
	ax.set_xlabel(r'Wavenumber (cm$^{-1}$)')
	ax.set_ylabel('%T data - lor fitted (single-peak fit)')
	plt.show()


def parse_args():
	
	#TODO: Make flag to toggle powerpoint and paper font settings
	parser = argparse.ArgumentParser()
	
	simulation_help = "Input file from transfer_matrix.py with transmittance, \
					   reflectance, absorptance data \
					   (not necessarily all three)."
	reference_help = "Input file from reference data (filmetrics, refractiveindex, etc.)"
	reflectance_help = "For testing: path for reflectance data downloaded from filmetrics."
	transmittance_help = "For testing: path for transmittance data downloaded from filmetrics"
	absorptance_help = "For testing: path for absorptance data downloaded from filmetrics"
	dispersion_help = "Plot dispersion curves from experimental angle-tuned data. \
					   Input path to dispersion data from polariton fitting."
	angle_help = "Plot cascading angle-tuned data on a single axis. \
						   Input path to angle-tuned data for a given sample concentration."
	splitting_help = "Results of nonlinear least squares to label dispersion plot. \
					  Input path to splitting data from polariton fitting."
	concentration_help = "Plots Rabi splitting versus concentration. Requires csv file \
						  generated by concentration_analysis module."
	residuals_help = "Plots residuals from spectra curve fitting. Needs residuals file."
	save_help = "Saves plot to a pdf instead of sending to the python viewer. \
				 Requires output path."
	units_help = "Sets the units for plotting. Valid units are: \
				  ev (eV), wn (wavenumber cm-1), wl (wavelength um)"


	parser.add_argument('-SIM', '--simulation', help=simulation_help)
	parser.add_argument('-TRN', '--transmission',help=transmittance_help)
	parser.add_argument('-RFL', '--reflection',help=reflectance_help)
	parser.add_argument('-ABS', '--absorbance', help=absorptance_help)
	parser.add_argument('-D', '--dispersion', help=dispersion_help)
	parser.add_argument('-T', '--angle', help=angle_help)
	parser.add_argument('--save_plot', help=save_help)
	parser.add_argument('-SP', '--splitting_results', help=splitting_help)
	parser.add_argument('-C', '--concentration', help=concentration_help)
	parser.add_argument('-R', '--residuals', help=residuals_help)
	parser.add_argument('units', type=str, nargs='*', help=units_help)
	
	return parser.parse_args()


def main():
	args = parse_args()
	if args.units:
		display_units = args.units[0].lower()

	if args.transmission:
		print("Using transmission reference data.")
		wl_T_data, T_data = reference_data(args.transmission)
	if args.reflection:
		print("Using reflection reference data.")
		wl_R_data, R_data = reference_data(args.reflection)
	if args.absorbance:
		print("Using absorbance reference data.")
		wl_A_data, A_data = reference_data(args.absorbance)


	if args.simulation:
		tmm_plots(args.simulation, args.save_plot)
# 	field_profile_data = sys.argv[2]

	if args.dispersion:
		#TODO: Use parameter strings to title plots
		dispersion_data = args.dispersion
		sample_name, params = pp.get_sample_params(dispersion_data)
		file_prefix = params[0] + '_' + params[1]
		if args.splitting_results:
			plot_dispersion(file_prefix, dispersion_data, display_units, args.splitting_results, args.save_plot)
		else:
			plot_dispersion(file_prefix, dispersion_data, display_units, args.save_plot)

	if args.angle:
		angle_data = args.angle
		sample_name, params = pp.get_sample_params(angle_data)
		
		file_prefix = params[0] + '_' + params[1]
		
# 		plot_spectra(file_prefix, angle_data, excitation=2171)
		plot_spectra(file_prefix, angle_data, save_dir=args.save_plot)
	
	if args.concentration:
		plot_splitting_concentration(args.concentration)
		
	if args.residuals:
		plot_curve_fit_residuals(args.residuals)
	
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


if __name__ == "__main__":
	main()

