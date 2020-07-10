import numpy as np
import scipy.fftpack as fft
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import polariton_processing as pp
import convert_unit
import csv


def gauss_model( omega, noise_=0., A=1, W=45, omega_0=2172):
	"""Gaussian line shape used for testing"""
	return A * np.exp(-((omega - omega_0) / W)**2) + noise_
	
def lorentz_model(omega, noise_=0., gamma=40., omega_0=2172):
	"""Lorentzian line shape used for testing"""
	return 1 / np.pi * (gamma**2 / ((omega - omega_0)**2 + gamma**2)) + noise_

def kramers_kronig(model_):
	"""Testing still"""

def write_refractive(frequency, real_n, imag_n, output_path):
	"""Write refractive index data to csv."""
	wavelength = convert_unit.wavenum_to_wavelen(frequency)
	
	file_name = output_path + 'test_hilbert_transform.csv'
	with open(file_name, 'w', newline='') as csvfile:
		csvwriter = csv.writer(csvfile)
		csvwriter.writerow(['"Wavelength, µm"', '"n"', '"k"'])
		for i, x in enumerate(wavelength):
			csvwriter.writerow([wavelength[i], real_n[i], imag_n[i]])
	return 0


def transform_absorbance(data_file, concentration, cavity_len, lbound_=2000, ubound_=2500):
	"""Perform Hilbert transform on FTIR absorbance data.
	   Return absorbance and transformed data."""

	omega_full, absorbance_full = pp.get_data(data_file)
	omega, absorbance = pp.truncate_data(omega_full, absorbance_full, lbound_, ubound_)
	extinction = absorbance / (concentration * cavity_len) # absorbance / (concentration * cavity length)
	transform = 1.33 - fft.hilbert(extinction)
	
	return omega, transform, extinction
	

def main():
	
	#Lower, upper bound on raw and fake data
	lbound = 500
	ubound = 4000
# 	n_points = 1000
	# First we need to get raw absorbance data.
	data_file = '/Users/garrek/projects/raw_data/absorbance/DPPA_in_DMF_2.0M_Abs_2.0M_DPPA_16scans_2.0cm.csv'
	output_path = '/Users/garrek/projects/pistachio/data/refractive_index_data/'
	hilbert_results = transform_absorbance(data_file, 2.0, 25., lbound, ubound)
	write_refractive(*hilbert_results, output_path)


# 	fig, ax = plt.subplots(3, sharex=True)
	fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
	
	omega, n_real, n_imag = hilbert_results
	ax1.plot(omega, n_real, color='darkblue')
	ax2.plot(omega, n_imag, color='maroon')

	ax1.spines['bottom'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax1.xaxis.tick_top()
	ax1.tick_params(labeltop=False)
	ax2.xaxis.tick_bottom()
	
	ax1.tick_params(axis='both', direction='in', right=True)
	ax2.tick_params(axis='both', direction='in', right=True)
	
	d = 0.015
	kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
	ax1.plot((-d, +d), (-d, +d), **kwargs)
	ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
	kwargs.update(transform=ax2.transAxes)
	ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
	ax2.plot((1-d, 1+d), (1-d, 1+d), **kwargs)
	

	ax2.set_xlabel(r'Wavenumber (cm$^{-1}$)', fontsize=14)
	fig.text(0.008, 0.5, 'Refractive Index', va='center', rotation='vertical', fontsize=14)
	
	save_file = '/Users/garrek/Documents/project_graphics/refractive_index_test.pdf'
	fig.savefig(save_file, bbox_inches='tight')
	
# 	ax.set_title('Raw Data and Hilbert transform')

	# A bunch o' testing with random noise and models etc.
	# CONCLUSION: Sci-py hilbert transform package might actually work.
# 
# 	noise = np.random.normal(0., 0.003, n_points)
# # 	noise = 0.
# 	omega_fake = np.linspace(lbound, ubound, n_points)
# 	gauss = gauss_model(omega_fake, noise, A=0.3, W=30, omega_0=2132) \
# 			+ gauss_model(omega_fake, noise, A=0.7, W=30, omega_0=2172) \
# 			+ gauss_model(omega_fake, noise, A=0.4, W=45, omega_0=2255)
# 	not_gauss = fft.hilbert(gauss)
# 	
# 	lorentz = lorentz_model(omega_fake, noise)
# 	not_lorentz = fft.hilbert(lorentz)
# 	
# 	ax[1].scatter(omega_fake, gauss, s=1)
# 	ax[1].scatter(omega_fake, not_gauss, s=1)
# 	ax[1].set_title('Testing gaussian line shape model')
# 	ax[2].scatter(omega_fake, lorentz, s=1)
# 	ax[2].scatter(omega_fake, not_lorentz, s=1)
# 	ax[2].set_title('Testing lorentz line shape model')
# 	
# 	plt.show()

if __name__ == "__main__":
	main()