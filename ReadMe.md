# Pistachio

Author: Garrek Stemo\
Date Created: July 5, 2019\
Updated: July 11, 2019\

Pistachio is a transfer matrix program based on Optical Waves in Layered Media by Pochi Yeh.
The goal is to create a versatile program that can handle an arbitrary number of different types of 
layers while still being easy to read and configure.


## Data directory

Data in the data directory are downloaded from refractiveindex.info unless stated otherwise.

## Config files

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
