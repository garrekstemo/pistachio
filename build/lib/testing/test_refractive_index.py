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
import kramers_kronig as kk
import transfer_matrix as tmm


def test_interpolation(yaml_filename):
	"""
	Test making new data points and interpolating that data in
	transfer_matrix module.
	"""
		
	layer_num = 1  # Change layer number to test different layer
	
	# Get default directory with a bunch of test files
	device = None
	with pkg_resources.path('default', yaml_filename) as yaml_data:
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


def test_transform(refractive_filename):
	"""
	Test Hilbert transform on imaginary part of refractive index data and compare
	to reference data.
	"""

	test_layer = None
	with pkg_resources.path('default', refractive_filename) as refractive_file:

		num_points = 100000
		min_wl = 2.0
		max_wl = 20.0
		test_layer = tmm.Layer('material', num_points, min_wl, max_wl, 0.0)
		test_layer.get_data_from_csv(refractive_file)
	
	original_wl = test_layer.wavelengths
	original_ri = test_layer.refractive_index
	original_ec = test_layer.extinction_coeff
	
	transform_ri = 1.329 - fft.hilbert(original_ec)
	
	fig, ax = plt.subplots()
	
	ax.plot(original_wl, original_ri, color='darkblue')
	ax.plot(original_wl, transform_ri)
	
	plt.show()
	

def main():
	
# 	test_interpolation('air_methanol10um_air.yaml')
# 	test_interpolation('fabry-perot_test.yaml')
# 	test_interpolation('fabry-perot_methanol_test.yaml')
# 	test_transform('Methanol.csv')

	
if __name__=='__main__':
	main()