import numpy as np
import scipy.fftpack as fft
from scipy import constants
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import polariton_processing as pp
import convert_unit
import csv



def square_root(x, x0, y0):
	"""Square root used for concentration-dependence.
	x0 is an optional normalization factor."""
	return [ y0 + np.sqrt(x_i / x0) for x_i in x]


def gaussian(w, amp, w_0, gamma):
	"""Gaussian function usually used for line shape"""
	return amp / (gamma * np.sqrt(2*np.pi)) * np.exp(- ((w - w_0)**2 / (2 * gamma**2)))


def lorentzian(w, amp, w_0, fwhm):
	"""Lorentzian function usually used for line shape"""
	return amp / np.pi * (fwhm/2)**2 / ((w - w_0)**2 + (fwhm/2)**2)


def pseudo_voigt(w, amp, w_0, gamma, m):
	"""
	Linear combination of Lorentzian and Gaussian distribution functions.

	w = frequency variable
	amp = Amplitude
	w_0 = peak center position
	gamma = phenomenological damping factor
	m = weight factor
	"""
	return m * gaussian(w, amp, w_0, gamma) + (1 - m) * lorentzian(w, amp, w_0, gamma)


def p(w, a, w_0, gamma):
	"""
	This perturbation function p(w) introduces an  asymmetry to Gaussian and Lorentzian line shapes.

	a=0 is the unperturbed (symmetric) system.
	a>0  corresponds to a tailing or skewing towards positive frequencies.
	a<0 corresponds to skewing in the direction of negative frequencies.

	Perform a substitution w -> w*p(w) in the Gaussian or Lorenztian functions
	to introduce this asymmetry:

	G(w * p(w)), L(w * p(w)) or 
	f(w, a) = m * G(w * p(w)) + (1-m) * L(w * p(w)).
	"""
	return 1 - a * (w - w_0) / gamma * np.exp(- (w - w_0)**2 / (2 * (2*gamma)**2))


def asym_voigt(w, amp, w_0, gamma, a, m):
	w = w * p(w, a, w_0, gamma)
	return pseudo_voigt(w, amp, w_0, gamma, m)


def double_asym_voigt(w, amp1, w1, gamma1, a1, m1, amp2, w2, gamma2, a2, m2):
	w1 = w1 * p(w1, a1, w1, gamma1)
	w2 = w2 * p(w2, a2, w2, gamma2)
	return asym_voigt(w, amp1, w1, gamma1, a1, m1) + asym_voigt(w, amp2, w2, gamma2, a2, m2)


def cavity_mode_energy(angle, E0, n_eff):
	return E0 / np.sqrt(1 - (np.sin(angle) /n_eff)**2)


def coupled_energies(theta, E0, Ee, V, n_eff, branch=0):
	"""Eigen-energies of coupled-state Hamiltonian."""
	Ec = E0 / np.sqrt(1 - (np.sin(theta) /n_eff)**2)
	if branch == 0:
		E_coupled = 0.5*(Ee + Ec) - 0.5*np.sqrt(4*V**2 + (Ee - Ec)**2)
	elif branch == 1:
		E_coupled = 0.5*(Ee + Ec) + 0.5*np.sqrt(4*V**2 + (Ee - Ec)**2)
	return E_coupled


def kramers_kronig(data_file, concentration, cavity_len, bounds=(-np.inf, np.inf), background=1.0):
	"""Rescale FTIR absorbance data and perform Hilbert transform.
	   Return transformed data."""

	omega_full, absorbance_full = pp.get_data(data_file)
	omega, absorbance = pp.truncate_data(omega_full, absorbance_full, bounds[0], bounds[1])
	extinction = absorbance / (concentration * cavity_len)
	transform = background - fft.hilbert(extinction)
	
	return omega, transform, extinction


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
	

# =============== Unit Conversions =============== #



def deg_to_rad(angles):
	"""Convert degrees to radians."""
	angles = [a * np.pi/180 for a in angles]
	return angles

def wavenumber_wavelength(wavenum):
	"""cm^-1 to micrometers"""
	wavelength = 10**4 / wavenum * 1.0
	return wavelength
	

# def wavelength_to_wavenumber(micrometers):
# 	"""micrometers to cm^-1"""
# 	wavenumbers = 10**4 / wave

def joules_to_ev(joules):
	ev = joules / constants.elementary_charge
	return ev

def wavenum_to_joules(wavenum):
	"""cm^-1 to photon energy"""
	cm_to_m = 1/100
	joules = constants.h * constants.c * (wavenum / cm_to_m)
	return joules

def wavenum_to_ev(wavenum):
	"""cm^-1 to eV units"""
	energy_J = wavenum_to_joule(wavenum)
	energy_ev = joule_to_ev(energy_J)
	return energy_ev

def ev_to_wavenum(energy_ev):
	"""Convert eV to cm^-1"""
	wavenumber = energy_ev / (constants.h * constants.c) * constants.elementary_charge / 100
	return wavenumber

def set_units(unit_data, current_units, set_units):
	#TODO: Handle arbitrary units from input data.
	# Assume input unit_data in cm^-1
	
# 	unit_data = np.array([unit_data])
	cm_to_m = 1/100
	energy = [wavenum_to_joules(d) for d in unit_data]
	energy_to_ev = [joules_to_ev(en) for en in energy]
	wavenumber_to_wavelength = [wavenumber_wavelength(d) for d in unit_data]

	if current_units == 'cm-1':
		if set_units == 'ev':
			new_units = energy_to_ev
		elif set_units == 'cm-1':
			new_units = unit_data
		elif set_units == 'um':
			new_units = wavenumber_wavelength

	elif current_units == 'ev':
		if set_units == 'ev':
			new_units = unit_data
		elif set_units == 'cm-1':
			new_units = [ev_to_wavenum(k) for k in unit_data]
		elif set_units == 'um':
			print("Why in the world would you want to convert to wavelength?")
			new_units = unit_data
	elif current_units == 'um':
		if set_units == 'ev':
			print('not right now')
		elif set_units == 'cm-1':
			new_units = wavenumber_wavelength(unit_data)

	return new_units


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