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

def write_refractive(frequency, real_n, imag_n, output_path, file_str):
	"""Write refractive index data to csv."""
	wavelength = convert_unit.wavenum_to_wavelen(frequency)
	
	file_name = output_path + file_str
	with open(file_name, 'w', newline='') as csvfile:
		csvwriter = csv.writer(csvfile)
		csvwriter.writerow(['"Wavelength, µm"', '"n"', '"k"'])
		for i, x in enumerate(wavelength):
			csvwriter.writerow([wavelength[i], real_n[i], imag_n[i]])
	print("Wrote real, imaginary refractive index data to", file_name)

def kramers_kronig(data_file, concentration, cavity_len, bounds=(-np.inf, np.inf), background=1.0):
	"""Rescale FTIR absorbance data and perform Hilbert transform.
	   Return transformed data."""

	omega_full, absorbance_full = pp.get_data(data_file)
	omega, absorbance = pp.truncate_data(omega_full, absorbance_full, bounds[0], bounds[1])
	extinction = absorbance / (concentration * cavity_len)
	transform = background - fft.hilbert(extinction)
	
	return omega, transform, extinction
	

def main():
	
	#Lower, upper bound on raw and fake data
	lbound = 500
	ubound = 4000
	L_cav = 9.3  # microns
# 	n_points = 1000

	output_path = '/Users/garrek/projects/pistachio/data/refractive_index_data/'

	# First we need to get raw absorbance data.
# 	data_file = '/Users/garrek/projects/raw_data/absorbance/DPPA_in_DMF/200715/Abs_DPPA_neat_BaF2_12umSpacer_16scans_2.0cm.csv'
# 	data_file = '/Users/garrek/projects/raw_data/absorbance/DPPA_in_DMF/200715/Abs_DPPA_DMF_0.285_BaF2_12umSpacer_16scans_2.0cm.csv'
# 	data_file = '/Users/garrek/projects/raw_data/absorbance/DPPA_in_DMF/200715/Abs_DPPA_DMF_0.913_BaF2_12umSpacer_16scans_2.0cm.csv'
	output_path = '/Users/garrek/projects/pistachio/data/refractive_index_data/'
# 	outfile_str = 'DPPA_DMF_0.285_BaF2_12umSpacer_L9.5.csv'
# 	hilbert_results = kramers_kronig(data_file, concentration=4.64, cavity_len=L_cav, bounds=(lbound, ubound), background=1.551)
# 	hilbert_results = kramers_kronig(data_file, concentration=0.285, cavity_len=L_cav, bounds=(lbound, ubound), background=1.487)
# 	hilbert_results = kramers_kronig(data_file, concentration=0.913, cavity_len=L_cav, bounds=(lbound, ubound), background=1.495)
	
	# 2.24M
# 	data_file = '/Users/garrek/projects/raw_data/absorbance/DPPA_in_DMF/200715/Abs_DPPA_DMF_2.240_BaF2_12umSpacer_16scans_2.0cm.csv'
# 	outfile_str = 'DPPA_DMF_2.24M_BaF2_12umSpacer_L9.5.csv'	
# 	hilbert_results = kramers_kronig(data_file, concentration=2.24, cavity_len=L_cav, bounds=(lbound, ubound), background=1.51)

	# 3.027M
# 	data_file = '/Users/garrek/projects/raw_data/absorbance/DPPA_in_DMF/200715/Abs_DPPA_DMF_3.027_BaF2_12umSpacer_16scans_2.0cm.csv'
# 	outfile_str = 'DPPA_DMF_3.027M_BaF2_12umSpacer_L9.5.csv'	
# 	hilbert_results = kramers_kronig(data_file, concentration=3.027, cavity_len=L_cav, bounds=(lbound, ubound), background=1.519)

	# 4.01M
	data_file = '/Users/garrek/projects/raw_data/absorbance/DPPA_in_THF/200803/Abs_DPPA_THF_4_CaF2_12umSpacer_16scans_2.0cm.csv'
	outfile_str = 'DPPA_THF_3.97M_CaF2_12umSpacer_L9.3.csv'	
	hilbert_results = kramers_kronig(data_file, concentration=3.97, cavity_len=L_cav, bounds=(lbound, ubound), background=1.506)

# 	write_refractive(*hilbert_results, output_path, outfile_str)
	


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
	
# 	save_file = '/Users/garrek/Documents/project_graphics/refractive_index_test.pdf'
# 	fig.savefig(save_file, bbox_inches='tight')
	
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
	plt.show()

if __name__ == "__main__":
	main()