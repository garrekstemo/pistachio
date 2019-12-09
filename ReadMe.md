# Pistachio

Author: Garrek Stemo\
Affiliation: Nara Institute of Science and Technology\
Date Created: July 5, 2019\
Updated: December 5, 2019

Pistachio is a suite of software for analysis of strongly coupled light-matter phenomenon. 
It includes a transfer matrix program based on Optical Waves in Layered Media by Pochi Yeh 
and Lorentzian fitting algorithms for generating dispersion curves, both from experimental data and simulated.


## Table of Contents

{{TOC}}


## Transfer Matrix Method

### Data Directory

The data directory contains some refractive index `csv` files to be used for simulations. They are taken from refractiveindex.info and filmetrics.com. The `data/out` directory is where simulation output data is stored, the user must input the file name and path as one of the command line arguments.

### Config files

Config files are stored in the .yaml format as a device with a configuration for each layer. At the top, the number of points and minimum and maximum wavelengths is specified, which are used to generate a list of wavelengths between the minimum and maximum for the simulation. Incident wave properties are specified next, including the incident angle, right-traveling wave amplitude and left-traveling wave amplitude. In the future, `theta_in` will be allowed to take a list of angles for angle-resolved simulations.

Each layer requires a material name, thickness (in nanometers), and a path to a file containing the refractive index for each wavelength. If there is no file, then the user must input the index of refraction and extinction coefficient; in this case, the user should put `None` next to `wavelength:`.

#### Example config file

```
num_points: 1000
min_wl: 1.
max_wl: 10.
wave:
	theta_in: 0.
	A0: 1
	B0: 0

layers:
	layer0:
	    material: SiO2
	    thickness: 0
	    param_path: "/../pistachio/data/layer0.csv"
    layer1:
	    material: Au
	    thickness: 10
	    param_path: "/../pistachio/data/layer1.csv"
    layer2:
	    material: Air
	    thickness: 1000
	    wavelength: None
	    index: 1.2
	    extinction: 0.0
	layer3
	    material: Au
	    thickness: 10
	    param_path: "/../pistachio/data/layer1.csv"
	layer4:
		material: SiO2
		thickness: 0
		param_path: "/../pistachio/data/layer0.csv"
```

Note that thickness must be given in nanometers. The param_path specifies
the path where refractiveindex.info data is stored, which must be saved as
a .csv file with wavelength, refractive index, and extinction coefficient
columns.


## Lorentzian fitting program

This is a bit of a mess right now, with way too many command-line flags and
code that is cobbled together. Some functions work and others need fine-tuning by the user. See the help (`polariton_processing.py -h`) documentation for run details.


## How to install Pistachio

Make sure you have python 3, numpy, scipy, and other dependencies installed.
The easiest way to do this is to install a package manager, like Anaconda.

Then go to the command line. Navigate to the directory where you want to put Pistachio, like `~/projects/`, using the `cd` command. Then enter

`git clone https://github.com/garrekds/pistachio .` 

## Running a transfer matrix simulation

This program runs from the command line. For example, if you want to run simulation with a p-wave field, you would navigate to the pistachio directory and enter

`python transfer_matrix.py -p config_files/file_name.yaml data/out/output_file.csv`

When in doubt, run `python transfer_matrix.py -h` to see the types and order of inputs.

## Plotting results

The output is a .csv file, so the user can use any plotting and analysis software. Basic plotting is provided via the included `plotting.py`. It will take the output file and generate a simple transmittance, reflectance, and absorbance plot. This may be made more sophisticated in the future.

## Things that don't work yet

The keen eye will notice a function that plots the field profile. This does not work yet; the field is discontinuous between boundaries when it should be continuous. I will fix it later.
