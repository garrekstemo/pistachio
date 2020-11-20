#! /anaconda3/bin/python

import csv
import os
import natsort as ns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns


def plot_spectra(dataframe, rawpath, ax, offset=0):
	"""
	Plot all spectra in a dataframe
	"""
	i = 0
	for index, row in dataframe.iterrows():
		row_path = str(row['Path'])
		if row_path != 'nan':
			datapath = rawpath + row_path
			data, absorb = get_angle_data_from_dir(datapath)
			for idx, spectrum in enumerate(data):
				if spectrum[0] == 0.0:
					angle, wavenumber, intensity = spectrum
					splitting = str(round(row['Splitting'], 1))
					concentration = str(round(row['Concentration'], 1))
					solvent = row['Solvent']
					label = r'$\Omega_R=$' + splitting + ', ' + solvent + ' C=' + concentration + 'M'
					ax.plot(wavenumber, intensity + i*offset, label=label)
		i+=1


def get_angle_data(simulation_path):
	"""
	Retrieve angles, transmission amplitudes, wavenumbers
	from transfer matrix calculations.
	"""
	angle_files = ns.realsorted(os.listdir(simulation_path))
	angle_data = []
	wavenumber_data = []
	transmission_data = []
	
	for i, angle in enumerate(angle_files):
		
		angle_file = os.path.join(simulation_path, angle)
		degree = float(angle[len('deg') : angle.find('_')])
		angle_data.append(degree)
		
		transmission = []
		wavenumber = []
		
		with open(angle_file, 'r', encoding='utf-8', newline='') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=',')
			header = next(csvreader)
			for row in csvreader:
				wavenumber.append(10**4 / float(row[0]))  # Convert from um to cm-1
				transmission.append(float(row[1]))
			wavenumber_data = wavenumber
			transmission_data.append(transmission)
	return [angle_data, wavenumber_data, transmission_data]


def more_points(data_array, num_points):
	
	multiple = int(num_points / len(data_array))
	mod = num_points % len(data_array)
	type(mod)
	
	new_array = []
	
	for a in data_array:
		new_array += multiple * [a]
	new_array += [data_array[-1]] * mod
	return new_array


def tmm_contour_plot(simulation_data, overlay=None):
	"""
	Produces a color plot using angle-resolved data to plot
	wavenumber vs angle with transmission as the z-axis.

	If save_plot is a filename and path, then the plot is saved in the format
	specified in the command line argument.
	"""

	angle, wavenumber, transmission = get_angle_data(simulation_data)
			
	angle = more_points(angle, len(wavenumber))
	transmission = more_points(transmission, len(wavenumber))
	Y, X = np.meshgrid(wavenumber, angle)
	fig, ax = plt.subplots()
	cbmin = 0.0
	cbmax = 0.1
	cbticks = np.linspace(cbmin, cbmax, 5)
	levels = np.linspace(cbmin, cbmax, 70)
	m = ax.contourf(X, Y, transmission, levels=levels)  # More levels = finer detail

	if overlay:
	# A csv file containing theta, upper/lower polariton wavenumber data from splitting fit.
		theta = []
		E_lp = []
		E_up = []
		with open(overlay, 'r', newline='') as f:
			csvreader = csv.reader(f, delimiter=',')
			for row in csvreader:
				theta.append(float(row[0]))
				E_lp.append(float(row[1]))
				E_up.append(float(row[2]))
		ax.plot(theta, E_lp, color='red', linestyle='dashed')
		ax.plot(theta, E_up, color='red', linestyle='dashed')
	
	# Plot formatting
	ax.set_ylim(1800, 2500)
	
	ax.tick_params(axis='both', labelsize=14)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
	ax.yaxis.set_minor_locator(ticker. AutoMinorLocator(5))
	ax.set_ylabel(r'Wavenumber (cm$^{-1}$)', fontsize=16)
	ax.set_xlabel('Angle (degrees)', fontsize=16)
# 	ax.set_title('4.0M DPPA in DMF')
	cbar = fig.colorbar(m, ticks=cbticks)
	cbar.set_label(label='Transmission (%)', size='large')
	cbar.ax.tick_params(labelsize='large')
	
# 	plt.show()
