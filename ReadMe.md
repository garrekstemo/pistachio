# Pistachio

Author: Garrek Stemo\
Affiliation: Nara Institute of Science and Technology\
Date Created: July 5, 2019\
Updated: November 20, 2020

Pistachio is a suite of software for analysis of coupled light-matter phenomenon.
It includes a transfer matrix program and fitting procedures for generating dispersion curves, both from experimental data and simulated.


## Table of Contents

{{TOC}}


## Transfer Matrix Method

### Introduction

The transfer matrix module computes transmission, reflection, and absorption for a multi-layer device using the transfer matrix method. Input waves can be s-wave or p-wave polarized. Mixed waves are not yet supported. The user creates a config file, which includes refractive index data for each layer in separate files. The program uses the multiprocess library for Python to parallelize transfer matrix calculations to distribute calculations for multiple incident angles to multiple processing cores.

This program has really only been tested thoroughly in Python 3.8.3. Your mileage may vary.

### Data Directory

The data directory contains some refractive index `csv` files to be used for simulations. They are taken from refractiveindex.info. To save transfer matrix results, the user must enter the output folder name as one of the command line arguments. The program will make a new folder based on parameters in a user-generated .yaml config file (example below).

### Config files

Config files are stored in the .yaml format with configurations for each material layer. At the top, the number of points, minimum, and maximum wavelengths are specified, which are used to generate a list of wavelengths between the minimum and maximum for the simulation. Incident wave properties are next, including the starting and ending incident angle (for angle-tuned measurements), right-traveling wave amplitude, and left-traveling wave amplitude. If `theta_i` (initial angle) and `theta_f` are different, then transfer matrix calculations are carried out for each angle in that range, with the number of angles determined by `num_angles`.

Each layer requires a material name, thickness (in nanometers), and a file containing the refractive index for each wavelength. This file must be placed in the `data/refractive_index_data` folder. If there is no file, then the user must input the index of refraction and extinction coefficient; in this case, the user should put `None` next to `wavelength` since these will be calculated based on max/min wavelengths and number of data points.

#### Example transfer matrix config file

```
num_points: 1000
min_wavelength: 1.0
max_wavelength: 10.0
wave:
	theta_i: 0.0
	theta_f: 30.0
	num_angles: 31
	A0: 1
	B0: 0

layers:
	layer0:
	    material: SiO2
	    thickness: 0
	    refractive_filename: "layer0.csv"
    layer1:
	    material: Au
	    thickness: 10
	    refractive_filename: "layer1.csv"
    layer2:
	    material: Air
	    thickness: 1000
	    wavelength: None
	    refractive_index: 1.0
	    extinction_coeff: 0.0
	layer3
	    material: Au
	    thickness: 10
	    refractive_filename: "layer1.csv"
	layer4:
		material: SiO2
		thickness: 0
		refractive_filename: "layer0.csv"
```

Note that thickness must be given in nanometers. The `refractive_filename` specifies the path where refractiveindex.info data is stored, which must be saved as a .csv file with wavelength, refractive index, and extinction coefficient columns.


### Running a transfer matrix simulation

This program runs from the command line. For example, if you want to run simulation with a p-wave field, you would navigate to the pistachio directory and enter

`python transfer_matrix.py -p config_files/file_name.yaml results`

When in doubt, run `python transfer_matrix.py -h` to see the types and order of inputs.


### How to name files and folders for experiments

In order to process data efficiently, it is important to have a consistent naming scheme. Parts of the program rely on file and directory naming conventions to batch process angle-resolved data. A directory containing angle-resolved .csv files shall have the naming convention

`concentration_solute_in_solvent`

For example:

`1.5M_WCO6_in_hexane`

The directory name is used to prevent the program from overwriting data when going through each angle.

Angle-resolved .csv files must contain the string `degNUM`where `NUM` is an integer. For example, `deg2` for incident angle 2 degrees. Inside this directory should also be an absorbance .csv file containing the target coupling band. This file should start with the string `Abs`.

Unfortunately, naming is done manually for most experimental setups so the user must take care to name their raw data carefully. No attempt is made by the program to guess misspellings, etc.

The program currently does not handle vacant cavity data, but this might be added in the future.


### Nonlinear least squares implementation

Upper and lower polariton experimental data are fitted to the positive and negative eigenvalues of the 2x2 Hamiltonian,

\\[
H = \left( \begin{matrix}
E_v & \hbar \Omega_i/2 \\
\hbar \Omega_i/2 & E_c
\end{matrix} \right),
\\]
where $\epsilon_v$ and $\epsilon_c$ are the vibration and cavity modes, respectively, and $\Omega_i$ is the vacuum Rabi splitting. Diagonalizing the Hamiltonian gives the eigenvalues corresponding the the upper and lower polariton modes[^2]. Further, we use the equation

\\[
E_c = E_0\left(1 - \frac{\sin^2(\theta)}{n_\text{eff}^2} \right)^{-1/2},
\\]

where $\theta$ is the light incident angle, $E_0$ is the cavity mode energy at 0 degree incidence, and $n_\text{eff}$ is the effective refractive index. Here we do not take into account the angle-dependence of the index of refraction. We use the SciPy least squares algorithm with initial guesses.


## Plotting and visualization

Experimental and simulated data are now processed in Jupyter notebooks included with this package. The Polarity Peak Analysis notebook is a workflow that takes the user through truncating spectra, determining a fitting model and parameters, batch fitting every angle for an angle-resolved experiment, and finally generating a dispersion curve, which is used to find the Rabi splitting parameter. A separate notebook uses uncoupled fringes and refractive index to determine cavity length. Detailed instructions are included in these notebooks.

## How to install Pistachio

Make sure you have Python 3, numpy, scipy, and other dependencies in the headers installed.
The easiest way to do this is to install a package manager, like Anaconda.

Then go to the command line. Navigate to the directory where you want to put Pistachio, like `~/projects/`, using the `cd` command. Then enter

`git clone https://github.com/garrekds/pistachio .` 

This is a Python package now, so you'll probably need to use something like `pip install --editable .` to install it. I don't exactly know how package installs work yet, but it works on my computer!

## References

[^1]: Yeh, P. Optical Waves in Layered Media. (Wiley, 2005).

[^2]: Skolnick, M. S., Fisher, T. A. & Whittaker, D. M. Strong coupling phenomena in quantum microcavity structures. Semicond. Sci. Technol. 13, 645–669 (1998).

[^3]: Weisstein, Eric W. "Lorentzian Function." From MathWorld--A Wolfram Web Resource. <http://mathworld.wolfram.com/LorentzianFunction.html>


