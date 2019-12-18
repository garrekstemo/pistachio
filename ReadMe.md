# Pistachio

Author: Garrek Stemo\
Affiliation: Nara Institute of Science and Technology\
Date Created: July 5, 2019\
Updated: December 12, 2019

Pistachio is a suite of software for analysis of strongly coupled light-matter phenomenon.
It includes a transfer matrix program based on Optical Waves in Layered Media by Pochi Yeh[^1] 
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
### Running a transfer matrix simulation

This program runs from the command line. For example, if you want to run simulation with a p-wave field, you would navigate to the pistachio directory and enter

`python transfer_matrix.py -p config_files/file_name.yaml data/out/output_file.csv`

When in doubt, run `python transfer_matrix.py -h` to see the types and order of inputs.

### Plotting results

The output is a .csv file, so the user can use any plotting and analysis software. Basic plotting is provided via the included `plotting.py`. It will take the output file and generate a simple transmittance, reflectance, and absorbance plot. This may be made more sophisticated in the future.

### Things that don't work yet

The keen eye will notice that there is a field profile function in the transfer matrix file. This does not work yet; the field is discontinuous at medium boundaries when it should be continuous. I'll fix this later.

## Lorentzian fitting and Dispersion Curves

This program processes and analyzes spectral data, specifically light-matter coupling FTIR experiments.

There are a few different flags that can be used to specify the type of analysis you want to do. They are: 

1. angle-resolved polariton + absorbance spectra -> dispersion curves (the main event)
2. bare cavity mode -> plot the spectrum and perform fitting
3. bare vibration mode -> plot the spectrum and perform fitting
4. single polariton spectrum -> plot and perform fitting

Command line arguments include, input file or directory, output directory, upper and lower wave number bound (for truncating data to perform fitting). An example input might look like

`python polariton_processing.py -B 2000 2400 -T /Users/../data/out/1.0M_WCO6_in_Hexane /Users/../data/out/` 


### How to name files and folders for experiments

In order to process data efficiently, it is important to have a consistent naming scheme. The `polariton_processing.py` program relies on file and directory naming consistency to batch process angle-resolved data. A directory containing angle-resolved csv files shall have the naming convention

`concentration_solvent_in_solute`

For example:

`1.5M_WCO6_in_hexane`

The directory name is used to prevent the program from overwriting data when going through each angle.

Angle-resolved csv files must contain the string `degNUM`where `NUM` is an integer. For example, `deg2` for incident angle 2 degrees. Inside this directory should also be an absorbance csv file containing the target coupling band. This file should start with the string `Abs`.

Unfortunately, naming is done manually for most experimental setups so the user must take care to name their raw data carefully. No attempt is made by the program to guess misspellings, etc.

The program currently does not handle vacant cavity data, but this might be added in the future.


### Angle-resolved polariton spectral data

Angle-resolved polariton data will include a Lorentzian[^3] double peak and multiple files representing each angle at which an experiment was performed. The directory may include absorbance data for the isolated target sample.

Each angle-resolved data must have the string `degx` somewhere in the file name, where `x` is the angle. The absorbance data should start with the string `Abs` so that can be properly processed.

The program will fit the absorbance data with a single-peak Lorentzian function and the angle-resolved data with a double-peak Lorentzian function. The peak wave numbers are then used to create dispersion curves for the upper and lower polaritons, as well as the vibration mode. The cavity mode is calculated using the 2x2 Hamiltonian

\\[
H = 
\left(
\begin{matrix}
\epsilon_v & \hbar \Omega_i/2 \\
\hbar \Omega_i/2 & \epsilon_c
\end{matrix}
\right),
\\]
where $\epsilon_v$ and $\epsilon_c$ are the vibration and cavity modes, respectively, and $\Omega_i$ is the vacuum Rabi splitting. Diagonalizing the Hamiltonian gives the eigenvalues corresponding the the upper and lower polariton modes that are found experimentally, leaving $\epsilon_c$ to be calculated[^2]. These calculations are still not quite perfectly implemented and are not to be trusted, but they do not produce Error messages!

More details on calculation methods later ....

Anyway, this option will output a csv will all the data needed to produce dispersion curves for the polariton modes, vibration mode, and cavity mode all in one plot. Proper plotting functionality coming later.

## Plotting

Plotting processed experimental data is relatively robust. Processed data from `polariton_processing.py` will have file names that can be read by `plots.py`, which can generate plots automatically. The user will need to manually adjust xlim and slim, as well as labels. This is a pain, but shouldn't be too hard to do.

## How to install Pistachio

Make sure you have python 3, numpy, scipy, and other dependencies installed.
The easiest way to do this is to install a package manager, like Anaconda.

Then go to the command line. Navigate to the directory where you want to put Pistachio, like `~/projects/`, using the `cd` command. Then enter

`git clone https://github.com/garrekds/pistachio .` 


# References

[^1]: Yeh, P. Optical Waves in Layered Media. (Wiley, 2005).

[^2]: Skolnick, M. S., Fisher, T. A. & Whittaker, D. M. Strong coupling phenomena in quantum microcavity structures. Semicond. Sci. Technol. 13, 645–669 (1998).

[^3]: Weisstein, Eric W. "Lorentzian Function." From MathWorld--A Wolfram Web Resource. <http://mathworld.wolfram.com/LorentzianFunction.html>


