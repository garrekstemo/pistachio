#!/anaconda3/bin/python

"""
Material properties and constants used in scripts.
Refractive index is valid for the wavelengths used in this research.
"""

class dppa(object):
	mol_mass = 275.204  # g / mol
	density = 1.277  # g / cm^3
	refractive = 1.551  # refractive index
	
class dmf(object):
	mol_mass = 73.1
	density = 0.944
	refractive = 1.43
	
class thf(object):
	mol_mass = 73.11
	density = 0.889
	refractive = 1.407
	
class toluene(object):
	mol_mass = 92.14
	density = 0.865
	refractive = 1.496

class phthalonitrile(object):
    mol_mass = 128.13  # g / mol
    density = 1.1  # g / cm^3
    refractive = 1.566  # WARNING: Estimate


class methanol(object):
    mol_mass = 32.042  # g / mol
    density = 0.791  # g / cm^3
    refractive = 1.329  # WARNING: Estimate
