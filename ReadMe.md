# Pistachio

Author: Garrek Stemo\
Affiliation: Nara Institute of Science and Technology\
Date Created: July 5, 2019\
Updated: July 8, 2020

Pistachio is a suite of software for analysis of coupled light-matter phenomenon.
It includes a transfer matrix program based on Optical Waves in Layered Media by Pochi Yeh[^1] 
and fitting algorithms for generating dispersion curves, both from experimental data and simulated.


## Table of Contents

{{TOC}}


## Transfer Matrix Method

### Introduction

The transfer matrix module computes transmission, reflection, and absorption for a multi-layer device using the transfer matrix method. Input waves can be s-wave or p-wave polarized. Mixed waves are not yet supported. The user creates a config file, which includes refractive index data for each layer in separate files. The program uses the multiprocess library for Python to parallelize transfer matrix calculations for angle-tuned calculations.

This program has really only been tested thoroughly in Python 3.8.3 on macOS 10.15.5. Your mileage may vary.

### Data Directory

The data directory contains some refractive index `csv` files to be used for simulations. They are taken from refractiveindex.info and filmetrics.com. To save transfer matrix results, the user must input the folder name as one of the command line arguments. The program will make a new folder based on parameters in a user-generated .yaml config file.

### Config files

Config files are stored in the .yaml format with configurations for each material layer. At the top, the number of points, minimum, and maximum wavelengths are specified, which are used to generate a list of wavelengths between the minimum and maximum for the simulation. Incident wave properties are next, including the starting and ending incident angle (for angle-tuned measurements), right-traveling wave amplitude, and left-traveling wave amplitude. If `theta_i` (initial angle) and `theta_f` are different, then transfer matrix calculations are carried out for each angle in that range, with the number of angles determined by `num_angles`.

Each layer requires a material name, thickness (in nanometers), and a file containing the refractive index for each wavelength. This file must be placed in the `data/refractive_index_data` folder. If there is no file, then the user must input the index of refraction and extinction coefficient; in this case, the user should put `None` next to `wavelength` since these will be calculated based on max/min wavelengths and number of datapoints.

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

### Plotting results

The output is a directory of .csv files for each angle (or just one if `theta_i` and `theta_f` are the same), so the user can use separate plotting and analysis software. Basic plotting is provided via the included `plots.py`. It will take the output file and generate simple transmittance, reflectance, and absorbance plots. This may be made more sophisticated analysis in the future.

### Things that don't work yet

The keen eye will notice that there is a field profile function in the transfer matrix file. This does not work yet. Further testing is needed.

There is a Kramers-Kronig function that kind of sort of works, but the user has to manually guide it a bit.

## Lorentzian fitting and Dispersion Curves

This program processes and analyzes spectral data, specifically light-matter coupling FTIR experiments. There will likely be significant changes to this in the future.

There are a few different flags that can be used to specify the type of analysis you want to do. They are: 

1. angle-resolved polariton + absorbance spectra -> dispersion curves (the main event)
2. bare cavity mode -> plot the spectrum and perform fitting
3. bare vibration mode -> plot the spectrum and perform fitting
4. single polariton spectrum -> plot and perform fitting

Command line arguments include, .yaml config file, input file or directory, and output directory, with the appropriate flags. An example input might look like the following:

`python polariton_processing.py -F /yaml_config.yaml -T /Users/../data/out/1.0M_WCO6_in_Hexane /Users/../data/out/` 

Here we have used the `-F` flag to denote a yaml config file containing the wavenumber bounds (for truncating data before doing a Lorentzian fit) and initial guesses for performing nonlinear least squares fitting on the dispersion data.


### How to name files and folders for experiments

In order to process data efficiently, it is important to have a consistent naming scheme. The `polariton_processing.py` program relies on file and directory naming consistency to batch process angle-resolved data. A directory containing angle-resolved .csv files shall have the naming convention

`concentration_solvent_in_solute`

For example:

`1.5M_WCO6_in_hexane`

The directory name is used to prevent the program from overwriting data when going through each angle.

Angle-resolved .csv files must contain the string `degNUM`where `NUM` is an integer. For example, `deg2` for incident angle 2 degrees. Inside this directory should also be an absorbance .csv file containing the target coupling band. This file should start with the string `Abs`.

Unfortunately, naming is done manually for most experimental setups so the user must take care to name their raw data carefully. No attempt is made by the program to guess misspellings, etc.

The program currently does not handle vacant cavity data, but this might be added in the future.

### Config file format

There is also a .yaml config file for analyzing polariton data. FTIR analysis is performed over a wide range of wavelengths, making any curve fitting impossible without truncating the raw data first. The user is spared from doing this themselves by simply including upper and lower bounds in the config file. The next item in the config file are initial guesses for the least squares fitting algorithm, which includes the 0-degree incident cavity mode energy, Rabi splitting parameter, refractive index, and vibrational mode energy. These energies can be put in units of eV or cm<sup>-1</sup> and are automatically converted for the fitting algorithm. The `units` parameter tells the program which units the initial guesses are in. These are then converted to cm<sup>-1</sup> because these are the standard units reported in FTIR machines, where these kinds of experiments are often performed. Yes, it's a little confusing that bounds are in inverse centimeters and the least squares guesses can be whatever you want, but this is a work in progress.

To simplify command-line arguments, the user can also specify an output directory here.

#### Example polariton config file

```
bounds:
	lower: 2000
	upper: 2400

least_squares_guesses:
	units: 'eV'
	E_cav_0: 0.275
	Rabi_splitting: 0.006
	refractive_index: 1.4
	E_exc: 0.275
	
data:
	output: '/Users/.../output_dir/'
```


### Dispersion curve generation

Angle-tuned polariton data will include a Lorentzian[^3] double peak and multiple files representing each angle at which an experiment was performed. The directory may include absorbance data for the isolated target sample.

Each angle-resolved data must have the string `degx` somewhere in the file name, where `x` is the angle as an integer. The absorbance data should start with the string `Abs` so that can be read properly.

The program will fit the absorbance data with a single-peak Lorentzian function and the angle-resolved data with a double-peak Lorentzian function. The peak wave numbers are then used to create dispersion curves for the upper and lower polaritons, as well as the vibration mode. 

### Nonlinear least squares implementation

Once the peak finding algorithm has been applied and the data reorganized into dispersion data, we fit the upper and lower polariton experimental data to the positive and negative eigenvalues of the 2x2 Hamiltonian,

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

where $\theta$ is the light incident angle, $E_0$ is the cavity mode energy at 0 degree incidence, and $n_\text{eff}$ is the effective refractive index. Here we do not take into account the angle-dependence of the index of refraction. We use the default SciPy least squares algorithm with initial guesses. For now, the output is in eV. In the future, the entire program may be made dimensionless to help with scaling and user unit preferences. The fitting gives reasonable results for these parameters, with slightly more variance in $n_\text{eff}$ with noisy and limited data.

Once fitting is complete, the polariton, cavity mode, and vibrational mode dispersion data is written to a file. The plotting program can use this output file to produce a ready-made plot for viewing.

## Plotting

Plotting processed experimental data is relatively robust. Processed data from `polariton_processing.py` will have file names that can be read by `plots.py`, which can generate some plots automatically. The user will need to manually adjust a few labels, positioning, and titles. This is a pain, but shouldn't be too hard to do. The process is slowly being automated to make pretty plots.

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


