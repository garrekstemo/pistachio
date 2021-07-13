#! /anaconda3/bin/python

import csv
import os
import natsort as ns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import seaborn as sns
import data_io
import pmath


def get_one_spectrum(data_series, raw_path):
	"""
	Plot spectrum in a pandas Series.
	"""
	label = ''
	wavenumber = np.array([])
	intensity = np.array([])
	data_set_path = str(data_series['Path'])
	
	if data_set_path != 'nan':
		data_path = raw_path + data_set_path
		dat, absorb = data_io.get_angle_data_from_dir(data_path)
		for idx, spectrum in enumerate(dat):
			if spectrum[0] == 0.0:
				angle, wavenumber, intensity = spectrum
				splitting = str(round(data_series['Splitting'], 1))
				concentration = str(round(data_series['Concentration'], 1))
				solvent = data_series['Solvent']
				label = r'$\Omega_R=$' + splitting + ', ' + solvent + ' C=' + concentration + 'M'
				
	return [wavenumber, intensity, label]
# 				ax.plot(wavenumber, intensity + offset, label=label)


def get_all_spectra(pd_data, raw_path):
	"""
	Plot all spectra in a dataframe
	"""
	
	plot_data = []
	if isinstance(pd_data, pd.DataFrame):
		for idx, row in pd_data.iterrows():
			# Each row in the pandas DataFrame is a Series.
# 			plot_offset = i * offset
			plot_info = get_one_spectrum(row, raw_path)
			plot_data.append(plot_info)
		
		return plot_data

	elif isinstance(pd_data, pd.Series):
		return [get_one_spectrum(pd_data, raw_path)]

	else:
		print('Not a valid pandas object.')



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


def plot_mode_dispersion(to_plot, ax=None, **kwargs):

    theta, energy = to_plot
    return ax.plot(theta, energy)
    

def save_plot_data(directory, plot_data, plot_name, col_names):
	"""
	Save experimental or fit data in easy-to-plot format.
	directory: should be the directory where raw data is stored (don't separate raw data from plot data).
	plot_data: the data you wish to plot later.
	plot_name: name of the save file.
	col_names: column names as a list of strings. Must match plot_data length.
	"""
	if ".csv" not in plot_name:
		plot_name = plot_name + ".csv"
	plot_file = directory + plot_name
	
	with open(plot_file, 'w', newline='') as f:
		
		csvwriter = csv.writer(f, delimiter=',')
		csvwriter.writerow(col_names)
		for row in plot_data:
			csvwriter.writerow(row)
		
		
		
		
		
		
	
