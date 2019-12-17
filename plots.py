#! /anaconda3/bin/python

import polariton_processing as pp
import argparse
import csv
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def TRA_plots(inputfile):
	wavelength = []
	transmittance = []
	reflectance = []
	absorptance = []
	field = []

	with open(inputfile, 'r') as f:
		reader = csv.reader(f)
		next(reader, None)
		for row in reader:
			wavelength.append(float(row[0]))
			transmittance.append(float(row[1]))
			reflectance.append(float(row[2]))
			absorptance.append(float(row[3]))
			
	print("Generating plots...")
	fig, axs = plt.subplots(3, 1, sharex=True)
	gs1 = gridspec.GridSpec(3, 1)
	gs1.update(wspace=0.025, hspace=0.005)

	ax = axs[0]
	ax.plot(wavelength, transmittance,
			color='b',
			linewidth=0.8,
			label="transfer matrix")
	if args.transmission:
		ax.plot(wl_T_data, T_data, linestyle="dashed", color='#FF5733', label="downloaded data")
	ax.set_ylabel('Transmittance %', fontsize=12)
	ax.tick_params(axis='both', labelsize=12)
# 	ax.set_xlim(1, 10)


	ax = axs[1]
	ax.plot(wavelength, reflectance,
			color='b',
			linewidth=0.8,
			label="transfer matrix")
	if args.reflection:
		ax.plot(wl_R_data, R_data, linestyle="dashed", color='#FF5733', label="downloaded data")
	ax.set_ylabel('Reflectance %', fontsize=12)


	ax = axs[2]
	ax.plot(wavelength, absorptance,
			color='b',
			linewidth=0.8,
			label="transfer matrix")
	if args.absorbance:
		ax.plot(wl_A_data, A_data, linestyle="dashed", color="#FF5733", label="downloaded data")
	ax.set_ylabel('Absorptance %', fontsize=12)
	ax.set_xlabel('Wavelength ($\mu$m)', fontsize=12)
	

	title = "Transfer Matrix Method"
	plt.suptitle(title, fontsize=18)
	plt.subplots_adjust(top=0.9)
# 	plt.tight_layout()

	if args.savedir:
		print("Saving figure as pdf")
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

def plot_spectra(file_prefix, spectra_file, excitation=None):
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
	i = len(angles) - 1
	while i >= 0:

		y_vals = intensities[i][1]
		y_vals = [y+y_offset for y in y_vals]

		deg_label = str(angles[i])
		ax.plot(wavenumbers, y_vals,
				color='black',
				linewidth=0.5,
				label=deg_label)
		y_offset += 0.25
		i-=1
		
	if excitation:
		ax.axvline(x=excitation)
	xy_pt1 = (2350, 0.09)
	xy_pt2 = (2350, 2.7)
	ax.annotate(r'$\theta$ = 0$\degree$', xy=xy_pt1)
	ax.annotate(r'$\theta$ = 20$\degree$', xy=xy_pt2)
	ax.set_xlim([2000, 2400])
	ax.set_ylim([0, 4])
	
	ax.set_xlabel(r'Wavenumber (cm$^{-1}$)')
	ax.set_ylabel('Transmission %')

	if args.savedir:
		file_name = file_prefix + '_' + 'cascade_plot.pdf'
		output_file = os.path.join(args.savedir, file_name)
		fig.savefig(output_file, bbox_inches='tight')
		print("Saved cascade plot to {}".format(output_file))
	else:
		plt.show()
	
	return 0
	

def plot_dispersion(file_prefix, dispersion_file):
	"""Takes dispersion curve file and plots wavenumber vs angle
	   for UP, LP, vibration mode, cavity mode"""
	   
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

	fig, ax = plt.subplots()
	
	vib_label = 'N=N=N  ' + str(int(vibration[0]))
	vib_xy = (angles[-1], vibration[0])
	
	cav_label = 'cavity dispersion'
	cav_xy = (angles[-1], cavity[-1])
	
	mark_size = 10
	
	ax.scatter(angles, up, s=mark_size)
	ax.scatter(angles, lp, s=mark_size)
	ax.plot(angles, vibration,
			linestyle='dashed',
			color='dimgray',
			label=vib_label)
	ax.plot(angles, cavity, linestyle='dashed', color='dimgray')
	
	ax.set_xticks(angles)
	
	ax.set_xlabel(r'Incident angle (deg)')
	ax.set_ylabel(r'Wavenumber (cm$^{-1}$)')
	
	ax.annotate(vib_label, xy=vib_xy, xytext=(-80, 3.), textcoords='offset points')
# 	ax.annotate(cav_label, xy=cav_xy, xytext=(-90, 0.), textcoords='offset points')
	
	if args.savedir:
		file_name = file_prefix + 'dispersion_curve.pdf'
		output_file = os.path.join(args.savedir, file_name)
		print(output_file)
		fig.savefig(output_file, bbox_inches='tight')
		print("Saved dispersion plot to {}".format(output_file))
	else:
		plt.show()
	
	return 0
	

def main():


	if args.simulation:
		TRA_plots(args.simulation)
# 	field_profile_data = sys.argv[2]

	if args.dispersion:
		#TODO: Handle Neat_DPPA type file names
		#TODO: Use parameter strings to title plots
		dispersion_data = args.dispersion
		sample_name, params = pp.get_sample_params(dispersion_data)
		file_prefix = params[0] + '_' + params[1]
		plot_dispersion(file_prefix, dispersion_data)

	if args.angles:
		angle_data = args.angles
		sample_name, params = pp.get_sample_params(angle_data)
		
		file_prefix = params[0] + '_' + params[1]
		
		plot_spectra(file_prefix, angle_data, excitation=2171)
# 		plot_spectra(file_prefix, angle_data)
	
	
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
	
	#TODO: Make flag to toggle powerpoint fonts and font sizes
	parser = argparse.ArgumentParser()
	
	simulation_help = "Input file from transfer_matrix.py with transmittance, \
					   reflectance, absorptance data."
	reference_help = "Input file from reference data (filmetrics, refractiveindex, etc.)"
	reflectance_help = "path for reflectance data downloaded from filmetrics"
	transmittance_help = "path for transmittance data downloaded from filmetrics"
	absorptance_help = "path for absorptance data downloaded from filmetrics"
	dispersion_help = "Use this to plot dispersion curves."
	angle_resolved_help = "Plot angle-resolved data on a single axis."
	save_help = "Saves plot to a pdf instead of sending to the python viewer.\
				 Requires output path."
	

	parser.add_argument('-SIM', '--simulation', help=simulation_help)
	parser.add_argument('-T', '--transmission',
						help=transmittance_help)
	parser.add_argument('-R', '--reflection',
						help=reflectance_help)
	parser.add_argument('-A', '--absorbance', help=absorptance_help)
	parser.add_argument('-D', '--dispersion', help=dispersion_help)
	parser.add_argument('-ANG', '--angles', help=angle_resolved_help)
	parser.add_argument('-S', '--savedir', help=save_help)
	
	args = parser.parse_args()

	if args.transmission:
		wl_T_data, T_data = reference_data(args.transmission)
	if args.reflection:
		wl_R_data, R_data = reference_data(args.reflection)
	if args.absorbance:
		wl_A_data, A_data = reference_data(args.absorbance)
	else:
		print("Using no reference data.\n")

	main()

