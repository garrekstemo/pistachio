import argparse
import csv
import os
import numpy as np
import polariton_processing as pp


def get_concentration_data(dir_str):
	"""Takes a directory path which should contain splitting
	   csv files from polariton_processing. Grabs the concentration from the
	   file name for each file, opens the file, grabs the Rabi splitting value.
	   Concentrations and Rabi splitting values are then put into an array and returned."""

	concentrations = []
	rabi_splittings = []
	solute = ''
	solvent = ''

	split_str = 'splitting_fit'
	for data_file in os.listdir(dir_str):
		if split_str in data_file:
			data_path = os.path.join(dir_str, data_file)
			sample_name, params = pp.get_sample_params(data_file)
			conc = float(params[0].strip('M'))
			solute = params[1]
			solvent = params[2]
			rabi = get_split_value(data_path)
			concentrations.append(conc)
			rabi_splittings.append(rabi)
	
	threaded = sorted(zip(concentrations, rabi_splittings))
	new_c, new_r = zip(*threaded)

	return solute, solvent, new_c, new_r

def get_split_value(splitting_file):
	"""Opens a file containing results of a splitting fit in polariton_processing.
	   Grabs the Rabi splitting value and returns it."""
	rabi_split = 0
	with open(splitting_file, 'r') as f:
		csvreader = csv.reader(f)
		for row in csvreader:
			if row[0].lower() == 'rabi':
				rabi_split = float(row[1])
	return rabi_split

def write_to_file(output_dir, solute, solvent, concentrations, splittings):
	"""Writes concentration and splitting data to a csv file in user-specified directory"""

	file_prefix = solute + '_in_' + solvent + '_concentration_dependence.csv'
	output = os.path.join(output_dir, file_prefix)
	with open(output, 'w') as f:
		filewriter = csv.writer(f, delimiter=',')
		header = ['Concentration (sqrt(M/M0))', 'Splitting Parameter (meV)']
		filewriter.writerow(header)

		i=0
		while i < len(concentrations):
			row = [concentrations[i], splittings[i]]
			filewriter.writerow(row)
			i+=1
	print("Wrote concentration dependence data to {}".format(output))
	return 0

def scale_concentration(concentration_list, c0=4.6):
	"""Performs math on raw concentration data"""
	
	new_list = []	
	for c in concentration_list:
		c = np.sqrt(c / c0)
		new_list.append(c)
		
	return new_list

def scale_splitting(splitting_list):
	"""Performs math on raw splitting data. 
	   USE ONLY IF SPLITTING IN EV."""
	
	new_list = []
	for i in splitting_list:
		i = i*1000
		new_list.append(i)
	return new_list
	
	
def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('conc_dir')
	parser.add_argument('out_file')
	return parser.parse_args()


def main():
	args = parse_args()

	solute, solvent, conc, splitting = get_concentration_data(args.conc_dir)
	conc = scale_concentration(conc)
# 	splitting = scale_splitting(splitting)
	write_to_file(args.out_file, solute, solvent, conc, splitting)


if __name__ == '__main__':
	main()