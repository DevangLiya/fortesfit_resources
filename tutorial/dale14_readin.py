# !!!!!!!!!!!!! ATTENTION: Edit line 23 before using this script !!!!!!!!!!!!!!!!
# This is the 09/06/2023 version.

import numpy as np
from scipy import integrate
from fortesfit.FortesFit_Settings import FortesFit_Cosmology as cosmo

""" FortesFit compliant readin module for the Dale et al. 2014 galaxy IR SED template library
"""

# ***********************************************************************************************
	
def	readtemplates(templateroot=None):
	""" Read in the templates that will be used by the associated readin routine
		
		The returned variable can have any user-defined form that is understood by the readin routine below.
	"""
	
	# Read Dale14 library from file
	if templateroot is None:
		# directory that works on my computer
		# Input the absolute path to the "Dale14" folder distributed with this tutorial
		templateroot = 

	if templateroot[-1] == '/':
		templateroot = templateroot[:-1]

	l_norm = 1e45

	AlphaD14       = np.linspace(0.06250, 4, 64)
	template_basic = np.loadtxt(f'{templateroot}/spectra.0.00AGN.dat')  #  Only the AGN-free template set is used here

	wavelength = template_basic[:,0]  #  base wavelengths are in microns
	seds = np.empty((len(wavelength),len(AlphaD14)))  #  Empty array which will store the scaled SEDs

	# To normalise, integrate the SED to get the total IR luminosity (8-1000 um)	
	intindex = (wavelength >= 8.0) & (wavelength <= 1000.0)

	for ialpha, alpha_val in enumerate(AlphaD14):
		
		template_orig = template_basic[:,ialpha+1]  # Single template corresponding to the alphad14 value

		# Dale 2014 templates are in arbitrary log nu*Fnu
		#  To integrate over wavelength, first convert to Flambda
		integrand = 10**(template_orig[intindex])/wavelength[intindex] # units are arbitrary, so simple division by rest wavelength
		lirtot = integrate.trapz(integrand,wavelength[intindex])
			
		seds[:,ialpha] = template_orig + (np.log10(l_norm) - np.log10(lirtot))  #  scale the template to log total IR luminosity of 45 erg/s
		
	templates = {'wavelength':wavelength, 'seds':seds, 'parameters':{'AlphaD14':AlphaD14, 'GalaxyIRLuminosity': np.log10(l_norm)}}
	
	return templates		
	

# ***********************************************************************************************

def		readin(parameters, redshift, templates=None):
	""" Given a specific parameter set and a redshift, return the Dale 2014 pure galaxy template in native units.

		The parameters are :
			GalaxyIRLuminosity: the 8-1000 um integrated luminosity from stellar systems in log erg/s/cm^2
			AlphaD14: The alpha shape parameter for the Dale 2014 library.

	"""
	
	# Templates are needed to process the readin for this class of model. Catch cases where the template is not supplied.
	if templates == None:
		raise ValueError('Templates must be supplied')
	# These templates are already normalised to L(8-1000) == 1e45 erg/s, in log nuLnu units.

	wave_basic = templates['wavelength']  #  base wavelengths are in microns
	wave = wave_basic*(1.0+redshift)  # Observed wavelength 

	param_basic = templates['parameters']['AlphaD14']
	# Identify the template index for the given value of alpha
	index = np.argmin(np.abs(param_basic - parameters['AlphaD14']))
	template_orig = 10**(templates['seds'][:,index])
	
	scale_factor = 10**(parameters['GalaxyIRLuminosity'] - templates['parameters']['GalaxyIRLuminosity'])
	
	lfactor = 4.0*np.pi*(cosmo.luminosity_distance(redshift).value*3.0856e24)**2.0 #area of spehere in cm^2
	observedflux = (template_orig*scale_factor/lfactor)/wave

	sed = {'observed_wavelength':wave,'observed_flux':observedflux}
	
	return sed

##############################################################################################################

############# testing #############
# lumin = 44
# alpha = 2.0
# redshift = 1.5

# templates = readtemplates()
# sed = readin({'AlphaD14': alpha, 'GalaxyIRLuminosity': lumin}, redshift, templates)
# waves = sed['observed_wavelength']
# flux = sed['observed_flux']

# print(f'Observed flux at {waves[261]/(1+redshift)} = {flux[261]}; templates output is {10**templates["seds"][:, 31][261]}')
# print(f'Observed flux at {waves[1261]/(1+redshift)} = {flux[1261]}; templates output is {10**templates["seds"][:, 31][1261]}')
# print(f'Observed flux at {waves[1464]/(1+redshift)} = {flux[1464]}; templates output is {10**templates["seds"][:, 31][1464]}')
# print(f'Observed flux at {waves[1471]/(1+redshift)} = {flux[1471]}; templates output is {10**templates["seds"][:, 31][1471]}')
###################################