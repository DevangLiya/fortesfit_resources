import numpy as np
from scipy import integrate # we will use this during normalisation
from fortesfit.FortesFit_Settings import FortesFit_Cosmology as cosmo # import cosmology from FortesFit

def readin(parameters, redshift, templates=None):
    """ Blackbody as a function of wavelength (um) and temperature (K).
 
    Returns units of lamdaFlambda
    Using actual wavelength function
    """
    # constants in micron units
    h = 6.62607015e-22
    c = 299792458000000.0
    k = 1.380649e-11
    T = parameters['Tbb'] # get temperature of the blackbody from user
    
    # we will normalise our SED to have a certain 8-1000 micron luminosity since it is
    # a useful physical quantity known to correlate with dust mass, SFR etc depending on the application
    user_norm = 10**parameters['logLIR']
    
    # define wavelength range in rest-frame
    # (I think) this samples a good range with enough granularity to catch turn-overs at typical temperatures 
    rest_waves = np.logspace(-2, 4, 1000)
    # also calculate observed frame wavelengths
    observed_waves = rest_waves*(1.0 + redshift)
    
    fl = (1/rest_waves**5)/(np.exp(h*c/(rest_waves*k*T)) - 1) # Flambda in arbitrary units
    lfl = rest_waves*fl # lFl is used for scaling etc
    
    # integrate to find the LIR
    intindex = (rest_waves >= 8.0) & (rest_waves <= 1000.0)
    lir = integrate.trapz(fl[intindex], rest_waves[intindex]) # !!!! NOTE: Integrate over Flambda !!!!
    I_norm = lfl*(user_norm/lir)
    
    # this luminosity will be diluted due to the distance. calculate that factor
    lfactor = 4.0*np.pi*(cosmo.luminosity_distance(redshift).value*3.0856e24)**2.0 #area of spehere in cm^2
    # apply the factor and convert to observed-frame flux
    observed_flux = (I_norm/lfactor)/observed_waves
    
    sed = {'observed_wavelength': observed_waves,'observed_flux': observed_flux}
    
    return sed