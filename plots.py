#! /anaconda3/bin/python

import argparse
import csv
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def TRA_plots(inputfile):
	wavelength = []
	transmittance = []
	reflectance = []
	absorptance = []
	field = []
	data = sys.argv[1]

	with open(data, 'r') as data:
		reader = csv.reader(data)
		next(reader, None)
		for row in reader:
			wavelength.append(float(row[0]))
			transmittance.append(float(row[1]))
			reflectance.append(float(row[2]))
			absorptance.append(float(row[3]))
			
	print("Generating plots...")
	fig, ax = plt.subplots(1, 1, sharex=True)
	gs1 = gridspec.GridSpec(3, 1)
	gs1.update(wspace=0.025, hspace=0.005)

# 	ax = axs[0]
	ax.plot(wavelength, transmittance, 'b-', label="transfer matrix")
	if args.trns:
		ax.plot(wl_T_data, T_data, linestyle="dashed", color='#FF5733', label="downloaded data")
	ax.set_ylabel('Transmittance %', fontsize=12)
	ax.tick_params(axis='both', labelsize=12)
	# 	ax.set_xlabel('wavelength ($\mu$m)', fontsize=20)
# 	ax.set_xlim(1, 10)
	# 	ax.legend()


# 	ax = axs[1]
# 	ax.plot(wavelength, reflectance, 'b-', label="transfer matrix")
# 	if args.refl:
# 		ax.plot(wl_R_data, R_data, linestyle="dashed", color='#FF5733', label="downloaded data")
# 	ax.set_ylabel('Reflectance %', fontsize=12)
# 
# 	ax = axs[2]
# 	ax.plot(wavelength, absorptance, 'b-', label="transfer matrix")
# 	if args.abso:
# 		ax.plot(wl_A_data, A_data, linestyle="dashed", color="#FF5733", label="downloaded data")
# 	ax.set_ylabel('Absorptance %', fontsize=12)
# 	ax.set_xlabel('wavelength ($\mu$m)', fontsize=12)
	

	title = "Etalon with 10 nm Au film, 10 micron air gap"
	# 	plt.suptitle(title, fontsize=24)
	plt.subplots_adjust(top=0.9)
	# 	plt.tight_layout()
	plt.show()


def reference_data(data_file):
	"""Gets reference data downloaded from
	websites. Filmetrics.com data are in nanometers"""
	wavelength = []
	Y = []
	unit = 1  # sets order of magnitude (nm or um)
	with open(data_file, 'r') as ref:
	
		reader = None
		if 'filmetrics' in str(ref):
			print('Filmetrics data. Units in nm')
			unit = 10**-3
			reader = csv.reader(ref, delimiter="\t")
		else:
			print("Refractiveindex.info data. Units in um")
			reader = csv.reader(ref)
		next(reader, None)  # Skip header
		for row in reader:
			wl = float(row[0]) * unit  # MAKE SURE UNITS ARE CORRECT
			wavelength.append(wl)
			Y.append(float(row[1]))
	return wavelength, Y


def main():
	TRA_data = sys.argv[1]
# 	field_profile_data = sys.argv[2]
	TRA_plots(TRA_data)
	
	
	# ===== Plotting FTIR data ===== #
# 	x_file = "/Users/garrek/Desktop/fpi0xwave.txt"
# 	y_file = "/Users/garrek/Desktop/trans0.txt"
# 	x_data = []
# 	y_data = []
# 	
# 	with open(x_file, 'r') as xf:
# 		lines = xf.readlines()
# 		for line in lines:
# 			x_data.append(float(line))
# 	with open(y_file, 'r') as yf:
# 		lines = yf.readlines()
# 		for line in lines:
# 			y_data.append(float(line))
# ============================== #


if __name__ == "__main__":
	
	#TODO: Make flag to toggle powerpoint fonts and font sizes
	parser = argparse.ArgumentParser()
	
	simulation_help = "Input file from transfer_matrix.py with transmittance, \
					   reflectance, absorptance data."
	reference_help = "Input file from reference data (filmetrics, refractiveindex, etc.)"
	reflectance_help = "path for reflectance data downloaded from filmetrics"
	transmittance_help = "path for transmittance data downloaded from filmetrics"
	absorptance_help = "path for absorptance data downloaded from filmetrics"
	
	parser.add_argument('sim_data', help=simulation_help)
	parser.add_argument('-T', '--trns',
						help=transmittance_help)
	parser.add_argument('-R', '--refl',
						help=reflectance_help)
	parser.add_argument('-A', '--abso', help=absorptance_help)
	
	args = parser.parse_args()

	if args.trns:
		wl_T_data, T_data = reference_data(args.trns)
	if args.refl:
		wl_R_data, R_data = reference_data(args.refl)
	if args.abso:
		wl_A_data, A_data = reference_data(args.abso)
	else:
		print("Using no reference data.\n")

	main()

