# Pistachio

Author: Garrek Stemo\
Affiliation: Nara Institute of Science and Technology\
Date Created: July 5, 2019\
Updated: November 28, 2019

Pistachio is a suite of software for analysis of strongly coupled light-matter phenomenon. 
It includes a transfer matrix program based on Optical Waves in Layered Media by Pochi Yeh 
and Lorentzian fitting algorithms for generating dispersion curves, both from experimental data and simulated.


## Table of Contents

{{TOC}}


## Transfer Matrix Method

### Data Directory

The data directory contains some refractive index `csv` files to be used for simulations. They are taken from refractiveindex.info and filmetrics.com. The `data/out` directory is where simulation output data is stored.

### Config files

Config files are stored in the .yaml format as a device with a configuration for each layer. Each layer requires a material name, thickness (in nanometers), and a path to a file containing the refractive index for each wavelength. If there is no file, then the user must input the index of refraction and extinction coefficient; in this case, the user should put `None` next to `wavelength:`.

```
num_points: 1000
min_wl: 0.2
max_wl: 10

bound1:
bound2:
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

Make sure you python 3, numpy, scipy, and other dependencies installed.
The easiest way to do this is to install a package manager, like Anaconda.

Then go to the command line. Navigate to the directory where you want to put Pistachio, like `~/projects/`, using the `cd` command. Then enter

`git clone https://github.com/garrekds/pistachio .` 

## Running a transfer matrix simulation

When in doubt, run `python transfer_matrix.py -h` to see the types and order of inputs.



