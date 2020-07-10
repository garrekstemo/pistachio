#!/usr/bin/env python
"""
Name: test_refractive_index
Description: Refractive index data downloaded from a database must 
			 be re-mapped to the datapoints (wavelengths) specified by the user
			 in the yaml config file. This remapping is done via the SciPy
			 interpolate function interp1d.
			 
			 For experiments, there is sometimes no refractive index information
			 available, so it must be calculated from absorbance data via a Hilbert
			 transform or Kramers-Kronig formula.
			 
			 This method tests:
			 
			 1. Interpolation remaps adequately to downloaded data, even if the user
			 	specifies max/min wavelengths that extends beyond the domain of the
			 	downloaded data (extrapolation).
			 
			 2. Tests that refractive index data calculated via transform matches the
			 	original downloaded data.
"""

import sys
import os
import importlib.resources as pkg_resources
import logging
import numpy as np
import scipy.fftpack as fft
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import kramers_kronig as kk
import transfer_matrix as tmm
import data.refractive_index_data


def test_interpolation(yaml_filename):
	"""
	Test making new data points and interpolating that data in
	transfer_matrix module.
	"""
		
	layer_num = 1  # Change layer number to test different layer
	
	# Get default directory with a bunch of test files
	device = None
	with pkg_resources.path(data.refractive_index_data, yaml_filename) as yaml_data:
		device = tmm.get_dict_from_yaml(yaml_data)

	layers = tmm.get_layers_from_yaml(device)
	single_layer = layers[layer_num]  
	
	original_ri = single_layer.refractive_index
	original_ec = single_layer.extinction_coeff
	original_wl = single_layer.wavelengths
	
	# Generate new data from yaml file
	min_wl = device['min_wavelength']
	max_wl = device['max_wavelength']
	num_points = device['num_points']
	wave = tmm.Wave(min_wl, max_wl, num_points)
	wave.make_wavelengths()
	new_wl = wave.wavelengths

	single_layer.make_new_data_points(new_wl)
	new_ri = single_layer.refractive_index
	new_ec = single_layer.extinction_coeff
	
	# Plot original and generated data
	fig, ax = plt.subplots()
	
	ax.plot(original_wl, original_ri, color='darkblue')
	ax.plot(original_wl, original_ec, color='maroon')
	ax.plot(new_wl, new_ri)
	ax.plot(new_wl, new_ec)
	
	plt.show()


def test_transform(refractive_filename, baseline=1.0):
	"""
	Test Hilbert transform on imaginary part of refractive index data and compare
	to reference data.
	"""

	test_layer = None
	num_points = 100000
	min_wl = 2.0
	max_wl = 20.0
	test_layer = tmm.Layer('material', num_points, min_wl, max_wl, 0.0)
	test_layer.get_data_from_csv(refractive_filename)
	material_name = refractive_filename.split('.')[0]
	
	original_wl = test_layer.wavelengths
	original_ri = test_layer.refractive_index
	original_ec = test_layer.extinction_coeff
	
	transform_ri = baseline - fft.hilbert(original_ec)
	
	fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

	my_fontsize = 18	
	
	ax1.set_title('Hilbert transform for {}'.format(material_name), fontsize=my_fontsize)
	ax1.plot(original_wl, original_ri, label="refractiveindex.info")
	ax1.plot(original_wl, transform_ri, label="transformed data")
	
	ax2.plot(original_wl, original_ec, color='black', label='extinction coefficient')
	
	ax1.spines['bottom'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax1.xaxis.tick_top()
	ax1.tick_params(labeltop=False)
	ax2.xaxis.tick_bottom()

	ax1.tick_params(axis='both', which='both', direction='in', right=True, top=True, labelsize=my_fontsize)
	ax2.tick_params(axis='both', which='both', direction='in', right=True, top=False, labelsize=my_fontsize)
	ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
	ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
	ax2.set_xlabel(r'Wavelength ($\mu$m)', fontsize=my_fontsize)
# 	ax2.set_ylabel('Real Refractive Index', fontsize=my_fontsize)

	# Make break in plot
	d = 0.015
	kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
	ax1.plot((-d, +d), (-d, +d), **kwargs)
	ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
	kwargs.update(transform=ax2.transAxes)
	ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
	ax2.plot((1-d, 1+d), (1-d, 1+d), **kwargs)


	fig.text(0.008, 0.5, 'Refractive Index',
			 va='center', 
			 rotation='vertical', 
			 fontsize=my_fontsize)
	text_str='baseline refractive index = {}'.format(baseline)
	ax1.annotate(text_str, (0.1, 0.82), xycoords='axes fraction',
				bbox=dict(boxstyle="round", fc='w'))
	ax1.legend()
	ax2.legend()
	
	plt.show()
	

def main():
	
# 	test_interpolation('air_methanol10um_air.yaml')
# 	test_interpolation('fabry-perot_test.yaml')
# 	test_interpolation('fabry-perot_methanol_test.yaml')
# 	test_transform('Methanol.csv', 1.329)
# 	test_transform('Ethanol.csv', 1.361)
# 	test_transform('DMMP.csv', 1.413)
	test_transform('Toluene.csv', 1.496)
if __name__=='__main__':
	main()