import numpy as np
from scipy import constants
import sys

def deg_to_rad(angles):
	"""Convert degrees to radians."""
	angles = [a * np.pi/180 for a in angles]
	return angles

def wavenum_to_wavelen(wavenum):
	"""cm^-1 to micrometers"""
	wavelength = 10000. / wavenum
	return wavelength
	
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
	
def set_units(unit_data, cur_units, set_units):
	#TODO: Handle arbitrary units from input data.
	# Assume input unit_data in cm^-1
	cm_to_m = 1/100
	energy = [wavenum_to_joules(d) for d in unit_data]
	energy_to_ev = [joules_to_ev(en) for en in energy]
	wavenumber_to_wavelength = [wavenum_to_wavelen(d) for d in unit_data]

	if cur_units == 'wn':
		if set_units == 'ev':
			new_units = energy_to_ev
		elif set_units == 'wn':
			new_units = unit_data
		elif set_units == 'wl':
			new_units = wavenumber_to_wavelength

	elif cur_units == 'ev':
		if set_units == 'ev':
			new_units = unit_data
		elif set_units == 'wn':
			new_units = [ev_to_wavenum(k) for k in unit_data]
		elif set_units == 'wl':
			print("Why in the world would you want to convert to wavelength?")
			new_units = unit_data
	elif cur_units == 'wl':
		print("Cannot convert from wavelength right now....")
		sys.exit()

	return new_units



def main():
	"""Nothing to execute"""
	
if __name__ == '__main__':
	main()