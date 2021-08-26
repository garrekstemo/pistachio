# Pistachio

Author: Garrek Stemo\
Affiliation: Nara Institute of Science and Technology\
Date Created: July 5, 2019\
Updated: July 27, 2021

Pistachio implements the Yeh 2x2 transfer matrix method [1].
As of July 2021, it has been almost completely refactored from a command line-only 
application to a Python package so that it can easily be
installed and classes and methods can be called from other programs.

Unit tests have been implemented to ensure the reproducibility of the code using Pytest. Unit tests and integration tests are still in active development with the goal of creating a program that stands up to scientific rigor and reproducibility with a test-driven approach that is lacking in the scientific programming community.

Pistachio uses the new type hinting features introduced in Python 3.9, which helps greatly with writing 
clean, reproducible code.
Type hinting also aids readability when distributing sometimes complex scientific code.

## Transfer Matrix Method

### Introduction

The transfer matrix module computes transmission and reflection for a multi-layered optical structure using the 2x2 transfer matrix method. The user may create a config file in the yaml format, which may include refractive index data for each layer in separate files. The program uses the multiprocess library for Python to parallelize transfer matrix calculations to distribute calculations for multiple incident angles to processing cores.

### Config files

Config files are stored in the yaml format with configurations for each material layer. At the top, the number of points, minimum, and maximum wavelengths are specified, which are used to generate a list of wavelengths between the minimum and maximum for the simulation. Incident wave properties include the starting and ending incident angle (for angle-tuned calculations), right-traveling wave amplitude, and left-traveling wave amplitude. If `theta_i` (initial angle) and `theta_f` are different, then transfer matrix calculations are carried out for each angle in that range, with the number of angles determined by `num_angles`.

Each layer requires a material name, thickness, and a file containing the refractive index for each wavelength. All values are assumed to be in SI units (i.e. meters instead of nanometers). The user can also specify the complex refractive index here or 
by directly accessing the attributes of the Layer class (demonstrated in the included
tutorials).

### Example transfer matrix yaml config file

```
num_points: 10000
min_wavelength: 1.6e-6
max_wavelength: 4.0e-6
wave:
    polarization: "s-wave"
    theta_i: 0.0
    theta_f: 30.0
    num_angles: 31
    A0: 1
    B0: 0
layers:
    layer0:
        material: CaF2
        thickness: 0.
        refractive_filename: "layer0.csv"
    layer1:
        material: Au
        thickness: 2.e-8
        refractive_filename: "layer1.csv"
    layer2:
        material: Air
        thickness: 1.e-5
        wavelength: 1.e-6
        refractive_index: 1.0003
        extinction_coeff: 0.0
    layer3:
        material: Au
        thickness: 2.e-8
        refractive_filename: "layer1.csv"
    layer4:
        material: CaF2
        thickness: 0.
        refractive_filename: "layer0.csv"
```

The `refractive_filename` specifies the path where refractive index data is stored, which must be saved as a .csv file with wavelength, refractive index, and extinction coefficient columns.


### Running a transfer matrix simulation

Pistachio may be run from the command line. For example, if you want to run simulation with a p-wave field, you would navigate to the pistachio directory and enter

`python transfer_matrix.py -p config_files/file_name.yaml results`

The `results` directory is any directory where you wish to write the results. When in doubt, run `python transfer_matrix.py -h` to see the types and order of inputs.

Classes and methods may also be called in other programs or in Jupyter Lab, which is now much more convenient and probably the preferred way.


## How to install Pistachio

Make sure you have Python 3, numpy, scipy, ruamel.yaml and other dependencies installed.
The easiest way to install dependencies is a package manager like [Anaconda](https://anaconda.org). Download the source code and put it in a folder of your choosing. This is a Python package, so you can use `pip install .` from the command line inside the top `pistachio` directory to install it once you download the source code.

Pistachio has now been published on [Pypi](https://pypi.org/project/pistachio-tm), so you can also install from there using `pip install pistachio-tm`.

You can test the code by simply navigating to the pistachio directory in the command line
and typing `pytest`, which will automatically find the test directory and files.


## References

[1] Yeh, P. Optical Waves in Layered Media. (Wiley, 2005)