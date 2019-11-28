# Pistachio

Author: Garrek Stemo\
Institution: Nara Institute of Science and Technology\
Date Created: July 5, 2019\
Updated: November 28, 2019

Pistachio is a suite of software for analysis of strongly coupled light-matter phenomenon. 
It includes a transfer matrix program based on Optical Waves in Layered Media by Pochi Yeh 
and Lorentzian fitting algorithms for generating dispersion curves, both
from experimental data and simulated.

## Transfer Matrix Method

### Data directory

Data in the data directory are downloaded from refractiveindex.info unless stated otherwise.

### Config files

Config files are stored in the .yaml format as a device with a specified
number of layers and configuration for each layer. For example, a device might
look like the following

num_layers: 2
layers:
    layer0:
        material: SiO2 
        thickness: 1000 
        param_path: "/Users/project/params/layer0.csv"
    layer1:
        material: Au
        thickness 35
        param_path: "/Users/project/params/layer1.csv"

Note that thickness must be given in nanometers. The param_path specifies
the path where refractiveindex.info data is stored, which must be saved as
a .csv file with wavelength, refractive index, and extinction coefficient
columns.


## Lorentzian fitting program

This is a bit of a mess right now, with way too many command-line flags and
code that is overall cobbled together. Some functions work and others need fine-tuning
by the user. See the help (-h) documentation for run details.
