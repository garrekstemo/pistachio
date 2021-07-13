import lmfit as lm
import numpy as np
import scipy.fftpack as fft
from scipy import constants
from scipy import optimize
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import csv


# =============== Fitting Procedures =============== #


def error_f_4param(x, theta, Elp_data, Eup_data):
	"""Here, we write the eigenvalues of the energy Hamiltonian for 
	   coupled cavity-vibration system. We minimize this with experimental data.
	"""
	
	E0, Ev, Rabi, n_eff = x
	
	Ec = cavity_mode_energy(theta, E0, n_eff)
	En = coupled_energies(Ec, Ev, Rabi, branch=0)
	Ep = coupled_energies(Ec, Ev, Rabi, branch=1)
	
	err1 = En - Elp_data
	err2 = Ep - Eup_data	
	err = np.concatenate((err1, err2))
			
	return err
	

def error_f_3param(x, theta, Ev, Elp_data, Eup_data):
	"""
	3-parameter error function for Rabi splitting fit.
	"""
	
	E0, Rabi, n_eff = x
	
	Ec = cavity_mode_energy(theta, E0, n_eff)
	En = coupled_energies(Ec, Ev, Rabi, branch=0)
	Ep = coupled_energies(Ec, Ev, Rabi, branch=1)
	
	err1 = En - Elp_data
	err2 = Ep - Eup_data	
	err = np.concatenate((err1, err2))
			
	return err

	
def error_f_2param(x, theta, Ev, n_eff, Elp_data, Eup_data):
	"""
	2-parameter error function for Rabi splitting fit.
	Here, we write the eigenvalues of the energy Hamiltonian for 
	coupled cavity-vibration system. We minimize this with experimental data.
	"""
	
	E0, Rabi = x
	Ec = cavity_mode_energy(theta, E0, n_eff)

	En = coupled_energies(Ec, Ev, Rabi, branch=0)
	Ep = coupled_energies(Ec, Ev, Rabi, branch=1)
	
	err1 = En - Elp_data
	err2 = Ep - Eup_data	
	err = np.concatenate((err1, err2))
			
	return err


def error_4peak(x, theta, Ev, polariton_data):
    """
    Simultaneous fitting of 4-peak strong coupling
    """
    E0, Rabi, n_eff = x
    Ec = cavity_mode_energy(theta, E0, n_eff)
    
    

def optimize_df(x, theta, E_up, E_lp):
	"""
	DEPRECATED: all fitting done in Jupyter notebooks.
	The derivative w.r.t. each variable for optimize_f function to
	   compute Jacobian in nonlinear least squares fit."""
	
	n = len(2*theta)
	m = len(x)
	jacobian = np.zeros((n, m))   
	
	E_0, E_vib, Rabi, n_eff = x
	E_cav = E_0 / np.sqrt(1 - np.sin(theta)**2 /n_eff**2)  #Can probably use the function
	
	dEc_dE0 = 1 / np.sqrt(1 - (np.sin(theta)**2 / n_eff**2))
	dEc_dn = E_0 * np.sin(theta)**2 / (n_eff*n_eff*n_eff * np.power((1 - (np.sin(theta)**2 / n_eff**2)), 1.5))
	root = np.sqrt(4*Rabi**2 + (E_vib - E_cav)**2)

	# Negative solution	
	dE_vib_neg = 0.5 - (0.5*(E_vib - E_cav) / root)
	dRabi_neg = -2*Rabi / root
	dE_0_neg = (0.5 + 0.5*(E_vib - E_cav) / root) * dEc_dE0
	dn_neg = (-0.5 - 0.5*(E_vib - E_cav) / root) * dEc_dn

	# Positive solution
	dE_vib_pos = 0.5 + (0.5*(E_vib - E_cav) / root)
	dRabi_pos = 2*Rabi / root
	dE_0_pos = (0.5 - 0.5*(E_vib - E_cav) / root) * dEc_dE0
	dn_pos = (-0.5 + 0.5*(E_vib - E_cav) / root) * dEc_dn

	dE_vib = np.concatenate((dE_vib_neg, dE_vib_pos))
	dRabi = np.concatenate((dRabi_neg, dRabi_pos))
	dE_0 = np.concatenate((dE_0_neg, dE_0_pos))
	dn_eff = np.concatenate((dn_neg, dn_pos))
	
	jacobian[:,0] = dE_0
	jacobian[:,1] = dRabi
	jacobian[:,2] = dn_eff
	jacobian[:,3] = dE_vib

	return jacobian


# =============== Unit Conversions =============== #


def deg_to_rad(angles):
	"""Convert degrees to radians."""
	angles = [a * np.pi/180 for a in angles]
	return angles

def wavenumber_wavelength(wavenum):
	"""cm^-1 to micrometers"""
	wavelength = 10**4 / wavenum * 1.0
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

def set_units(unit_data, current_units, set_units):
	#TODO: Vectorize this math
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


# =============== Fitting Models =============== #


def refractive_line(x, n1, n2, c1, c2):
	"""
	Calculates refractive index for a particular concentration (mol / L) x
	for neat liquid refractive indices n1 and n2 with concentrations c1 and c2.
	Index 2 is the solute and index 1 is the solvent, thus c1 is zero and
	c2 is the solute concentration at room temperature.
	Linearity can be assumed based experimental data obtained for this research
	and also from the following reference: (TODO: Include that paper I found.)
	"""
	return np.abs((n2 - n1) / (c2 - c1)) * x + n1


def square_root(x, x0, y0):
	"""Square root used for concentration-dependence.
	x0 is an optional normalization factor."""
	return [ y0 + np.sqrt(x_i / x0) for x_i in x]


def gaussian(w, amp, w_0, gamma):
	"""Gaussian function usually used for line shape"""
	return amp / (gamma * np.sqrt(2*np.pi)) * np.exp(- ((w - w_0)**2 / (2 * gamma**2)))


def lorentzian(w, amp, w_0, fwhm):
	"""Lorentzian function usually used for line shape"""
	return (amp / (2*np.pi)) * fwhm / ((w - w_0)**2 + (fwhm/2)**2)
	

def lorentz_imag(w, amp, w_0, gamma):
	return amp * gamma * w / ((w**2 - w_0**2)**2 + (w*gamma)**2)

	
def lorentz_real(w, amp, w_0, gamma, offset):
	return offset + amp * (w**2 - w_0**2) / ((w_0**2 - w**2)**2 + (w*gamma)**2)


def complex_lorentzian(w, amp, w_0, gamma):
	return amp / (w**2 - w_0**2 + 1j*w*gamma)


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
	See this article for details: https://arxiv.org/abs/1804.06083

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
	"""
	This function is somewhat redundant because you can just add (+) models 
	in lmfit according to the CompositeModel class. But feel free to use this function
	for testing or something.
	"""
	w1 = w1 * p(w1, a1, w1, gamma1)
	w2 = w2 * p(w2, a2, w2, gamma2)
	return asym_voigt(w, amp1, w1, gamma1, a1, m1) + asym_voigt(w, amp2, w2, gamma2, a2, m2)


# def asym_lorentzian()

def make_one_model(num, peaks, amplitude, sigma):
    """
    Build one Lorentzian model.
    """
    pref = 'l{}_'.format(num)
    mod = lm.models.LorentzianModel(prefix=pref)
    params = mod.make_params()
    params[pref+'center'].set(peaks[num])
    params[pref+'amplitude'].set(amplitude[num])
    params[pref+'sigma'].set(sigma[num])

    return mod, params


def make_model(peaks, amplitude, sigma):
    """
    Build composite model from make_one_model function.
    """

    model = None
    params = None
    for i in range(len(peaks)):
        add_model, add_params = make_one_model(i, peaks, amplitude, sigma)
        if model is None:
            model = add_model
            params = model.make_params()
        else:
            model = model + add_model
        params.update(add_params)

    return model, params


# =============== Experiment-Specific Calculations =============== #

def kramers_kronig(data_file, concentration, cavity_len, bounds=(-np.inf, np.inf), background=1.0):
	"""Rescale FTIR absorbance data and perform Hilbert transform.
	   Return transformed data."""

	absorbance = data_file
	extinction = absorbance / (concentration * cavity_len)
	transform = background - fft.hilbert(extinction)
	
	return transform


def write_refractive(frequency, real_n, imag_n, output_path, file_str):
	"""Write refractive index data to csv."""
	wavelength = 10**4 / frequency
	
	file_name = output_path + file_str
	with open(file_name, 'w', newline='') as csvfile:
		csvwriter = csv.writer(csvfile)
		csvwriter.writerow(['"Wavelength, µm"', '"n"', '"k"'])
		for i, x in enumerate(wavelength):
			csvwriter.writerow([wavelength[i], real_n[i], imag_n[i]])
	print("Wrote real, imaginary refractive index data to", file_name)
	

def cavity_mode_energy(theta, E0, n_eff):
    """
    Recall that E0 can be either a quantity on its own or
    E0 = 1 / n_eff * (mode / L_cav) / 100
    """
    return E0 / np.sqrt(1 - (np.sin(theta) /n_eff)**2)


def coupled_energies(Ec, Ev, Rabi, branch=0):
	"""Eigen-energies of coupled-state Hamiltonian."""

	if branch == 0:
		E_coupled = 0.5*(Ev + Ec) - 0.5*np.sqrt(Rabi**2 + (Ev - Ec)**2)
	elif branch == 1:
		E_coupled = 0.5*(Ev + Ec) + 0.5*np.sqrt(Rabi**2 + (Ev - Ec)**2)

	return E_coupled
	
	
def ultrastrong_RWA(Ec, Ev, Rabi, branch=0):
    """Eigenenergies for ultrastrong coupling applying the rotating wave approximation (RWA)"""
    
    root = np.sqrt(4*Rabi**2 + (Ec - Ev)**2)

    if branch == 0:
        E_coupled = Ec + Ev + root
    elif branch == 1:
        E_coupled = Ec + Ev - root
    elif branch == 2:
        E_coupled = -Ec - Ev + root
    elif branch == 3:
        E_coupled = -Ec - Ev - root
        
    return 0.5 * E_coupled
        

def ultrastrong_full(Ec, Ev, Rabi, branch=0):
    """Eigenenergies for full Hamiltonian for ultrastrong coupling"""

    D = Rabi**2 / Ev
    first = Ec**2 + Ev**2 + 4*D*Ev
    inner_root = np.sqrt(16*Ec*Ev*Rabi**2 + (8*D*Ev)**2 + 8*D*Ev**2 - 8 * D * Ec**2 * Ev + (Ec**2 - Ev**2)**2)
    
    if branch == 0:
        E_coupled = np.sqrt(first + inner_root)
    elif branch == 1:
        E_coupled = np.sqrt(first - inner_root)
    elif branch == 2:
        E_coupled = -np.sqrt(first - inner_root)
    elif branch == 3:
        E_coupled = -np.sqrt(first + inner_root)
    
    return E_coupled / np.sqrt(2)


def main():
    sns.set_theme(palette='dark')

    fig, ax = plt.subplots(2, 2, sharey=True, sharex=False)
    
    RWA_color = 'b'
    full_color = 'r'
    
    L_cav = 25 * 10**-6
    modes = np.array([1, 2, 3, 4])
    theta = np.linspace(0, np.pi/3, 100)
    theta_deg = theta * 180/np.pi
    length = np.linspace(4, 13, 100) * 10**-6
    length_um = length * 10**6
    Ev = 2171
    Rabi = 88
    n_eff = 1.53

    for m in modes:
        E0_var = m / (n_eff * length) / 100
        E0_fixed = m / (n_eff * L_cav) / 100
        Ec = cavity_mode_energy(0.0, E0_var, n_eff)
        Ec_th = cavity_mode_energy(theta, E0_fixed, n_eff)
        
        branch = 0
        while branch <= 1:
            RWA = ultrastrong_RWA(Ec, Ev, Rabi, branch)
            full = ultrastrong_full(Ec, Ev, Rabi, branch)
            
            RWA_th = ultrastrong_RWA(Ec_th, Ev, Rabi, branch)
            full_th = ultrastrong_full(Ec_th, Ev, Rabi, branch)
        
            ax[0, 0].plot(length_um, RWA, label='RWA', color=RWA_color)
            ax[0, 1].plot(length_um, full, label='Full Hamiltonian', color=full_color)
            
            ax[1, 0].plot(theta_deg, RWA_th, label='RWA', color=RWA_color)
            ax[1, 1].plot(theta_deg, full_th, label='Full Hamiltonian', color=full_color)
#             ax[1, 1].plot(theta_deg, full_th, label='Full Hamiltonian', color=full_color)
            branch += 1


#     ax[0].set_ylim(1700, 2400)

    ax[0, 0].set_ylabel(r'Wavenumber (cm$^{-1}$)')
    ax[1, 0].set_ylabel(r'Wavenumber (cm$^{-1}$)')
    ax[0, 0].set_xlabel(r'Thickness $\mu$m (deg)')
    ax[0, 1].set_xlabel(r'Thickness $\mu$m (deg)')
    ax[1, 0].set_xlabel(r'$\theta$ (deg)')
    ax[1, 1].set_xlabel(r'$\theta$ (deg)')
    ax[0, 0].set_title('RWA')
    ax[0, 1].set_title('Full Hamiltonian')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
	main()