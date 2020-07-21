import ares
import numpy as np
import matplotlib.pyplot as pl
import os
import scipy
from ares.physics.Constants import c, h, k_B, g_per_msun

path = os.path.realpath(os.path.dirname(__file__))
os.chdir(path)

def convert_SED(light_source):
    """
    (SynthesisModel object) -> 1 or 2 dimensional ndarray

    Converts SED data from ergs / s / angstrom / [depends] to ergs / s / Hz / [depends] and 
    returns the data array. The input must be a proper SynthesisModel object
    containing the 'data' and 'dwdn' attributes.
    """
    return light_source.data * (np.outer(light_source.dwdn, np.ones(light_source.data.shape[1])))


# TODO Adapt this to accomodate star burst or continuous formation. For now this is only continuous.
def get_L_nu(light_source, SFR):
    """
    (SynthesisModel object, number) -> 1 or 2 dimensional ndarray

    Calculates specific luminosity in ergs / s / Hz from source SED data 
    and returns the data array. The input must be a proper SynthesisModel object
    containing the 'data' and 'dwdn' attributes.

    SFR must be given in Msun / year.
    """
    temp = convert_SED(light_source)
    return temp * SFR







# Temporary tests
if __name__ == "__main__":
    src = ares.sources.SynthesisModel(source_sed='eldridge2009', source_Z=0.02)
    new_SED = convert_SED(src)

    pl.loglog(src.wavelengths, src.data[:,0], label = 'wavelength')
    pl.loglog(src.frequencies, new_SED[:,0], label = 'frequency')
    pl.legend()
    pl.show()

    L_nu = get_L_nu(src, 1)
    pl.loglog(src.frequencies, L_nu[:,0])
    pl.show()