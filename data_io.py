#! /anaconda3/bin/python

import os
import re
import sys
import csv
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import optimize
from scipy.interpolate import interp1d
from scipy import constants
from scipy import signal
import matplotlib.pyplot as plt
from ruamel_yaml import YAML
import pdb
import pmath

yaml = YAML()


# ========== Get paramaters and data from user inputs and files ========== #

def get_output_path_from_yaml(yaml_config):
	"""Gets data output path from yaml config file.
	   This saves the user a command-line argument.
	   You're welcome.'"""	
	   
	with open(yaml_config, 'r') as yml:
		config = yaml.load(yml)
	return config['data']['output']
	

def get_FTIR_data(spectral_data, data_format=None):
	"""
	Retrieve csv spectral data from JASCO FTIR output and store in array.
	"""

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
	NOT TESTED.
	"""

	wavenumber_list = np.array([])
	intensity_list = np.array([])
	with open(spectral_data, 'r', errors='ignore') as spectrum:
		csvreader = csv.reader(spectrum)
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


def from_tmm(directories, angle_index=0):
	"""
	Import several transfer matrix sim directories,
	pull out the specified angle file for each sim,
	load them up into a pandas datafrom, and
	export to csv for comparison later.
	
	directories must be in the form of a dictionary with
	headers 'Sample' and 'Path'.
	
	Used in TMM Explorer Jupyter notebook
	"""

	df = pd.DataFrame(data=directories)
	df['Path']
	out_df = pd.DataFrame()
	
	for i, p in enumerate(df['Path']):
		simlist = os.listdir(p)
		spectrum = os.path.join(p, simlist[angle_index])
		transmittance = pd.read_csv(spectrum, usecols=['Transmittance'])
		out_df[df['Sample'][i]] = transmittance['Transmittance']
		
		if i == 0:
			wavelength = pd.read_csv(spectrum, usecols=['Wavelength'])
			out_df['Wavenumber'] = 10**4 / wavelength['Wavelength']

	headers = ['Wavenumber'] + df['Sample'].tolist()
	out_df = out_df.reindex(columns=headers)
	out_df = out_df.sort_values(by='Wavenumber')
	return out_df


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



def main():
	"None"


if __name__ == '__main__':
	main()
