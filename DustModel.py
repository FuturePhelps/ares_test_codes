# Author: Felix Bilodeau-Chagnon

# Created on 2020/06/25 at 11:31 AM EST

# Based on the paper by Imara et al. (2018):
# "A Model Connecting Galaxy Masses, Star 
# Formation Rates, and Dust Temperatures
# across Cosmic Time"

# Purpose: Create a model of infrared
# emissions due to dust in galaxies

# Main application: Obtain an emission
# spectrum for dust given parameters
# from a galaxy simulation

import numpy as np
import scipy
from ares.physics.Constants import c, h, k_B, g_per_msun
from ares.util.ParameterFile import ParameterFile
from ares.physics import Cosmology

F_STAR = 0.25

# TODO Change the inputs because they don't match what ARES gives us
# Inputs needed: L_nu (3d array for galaxy population? maybe 2d array and make another class)
# frequencies, R_halo, M_gas, z (1darray)

# TODO I think I just pass in a GalaxyEnsemble object into this...
# This will make it easier to access everything I need in this class,
# and it makes more sense because it is the only one which simulates
# dusts mass etc. Probably means that we need to put this file in
# the populations directory. Then we can also use the SPS frequencies
# as defaults?

# TODO this needs a lot of rewriting now...
# Pass in the ParameterBundle from a GalaxyEnsemble instance? Or maybe just the histories? Not sure.

class DustModel:
    @staticmethod
    def validate_data(L_nu, frequencies, R_halo, M_dust, z):
        """
        Returns a boolean.

        Validates if the data provided is correct. If not, raises a ValueError.
        frequencies must correspond to the first dimension of L_nu, and all
        other data must correspond to the second dimension of L_nu.
        """
        if frequencies.shape[0] != L_nu.shape[0]:
            raise ValueError("frequencies and L_nu have incompatible shapes")

        if (isinstance(R_halo, np.ndarray)) and (R_halo.shape != L_nu.shape[1:]):
            raise ValueError("R_halo and L_nu have incompatible shapes")

        if (isinstance(M_dust, np.ndarray)) and (M_dust.shape != L_nu.shape[1:]):
            raise ValueError("M_dust and L_nu have incompatible shapes")

        if (isinstance(z, np.ndarray)) and (z.shape != L_nu.shape[1:]):
            raise ValueError("z and L_nu have incompatible shapes")

        return True

    def must_recalculate(self, attribute):
        """
        (str) -> boolean

        Determines whether a certain attribute must be recalculated (due to data
        updating).
        """
        return (not hasattr(self, attribute)) or self._modified[attribute]

    # Parameters needed: L_nu (ergs / s / Hz), frequencies (Hz), R_halo (cm), M_dust (solar masses), z 
    def __init__(self, L_nu, frequencies, R_halo, M_dust, z, cosm = None, **kwargs):
        """
        (2darray, 1darray, [flex], [flex], [flex], [flex], Cosmology, **kwargs) -> DustModel

        [flex] = 1darray or number

        Initializes a DustModel object, which is used to
        calculate dust emissions from existing galaxy simulation
        models.

        L_nu: specific luminosity at all frequencies provided [ergs / s / Hz]
        frequencies: in [Hz]
        R_halo: radius of dark matter halo [cm]
        M_dust: dust mass of the galaxy(ies) [solar masses]
        z: redshift(s)

        The data provided must be in units of ergs / s / Hz (power per frequency).

        All data and frequencies must be provided in ascending order of
        frequencies. The second dimension of the L_nu array can be any
        other property (such as different times, different galaxies, etc.)
        """

        # Connect with parameter file and cosmology is needed later
        # TODO Currently not in use
        self.pf = ParameterFile(**kwargs)
        self._cosm_ = cosm

        # Check that sizes match
        DustModel.validate_data(L_nu, frequencies, R_halo, M_dust, z)

        # Save the data and associated frequencies
        self._L_nu = L_nu
        self._frequencies = frequencies
        self._R_halo = R_halo
        self._M_dust = M_dust
        self._z = z

        # Used to determine if we need to recalculate properties
        self._modified = {}

    @property
    def L_nu(self):
        """
        (void) -> 2darray

        Returns the L_nu array currently stored in this instance.
        """
        if hasattr(self, '_L_nu'):
            return self._L_nu
        return None

    @property
    def frequencies(self):
        """
        (void) -> 1darray

        Returns the frequencies associated with the data currently
        stored in this instance.
        """
        if hasattr(self, '_frequencies'):
            return self._frequencies
        return None
    
    @property
    def R_halo(self):
        """
        (void) -> 1darray / num

        Returns the halo radius / radii currently stored in this
        instance.
        """
        if hasattr(self, '_R_halo'):
            return self._R_halo
        return None

    @property
    def M_dust(self):
        """
        (void) -> number / 1darray

        Returns the dust mass(es) currently stored in this
        instance.
        """
        if hasattr(self, '_M_dust'):
            return self._M_dust
        return None

    @property
    def z(self):
        """
        (void) -> number / 1darray

        Returns the redshift(s) currently stored in this instance.
        """
        if hasattr(self, '_z'):
            return self._z
        return None

    def set_data(self, L_nu, frequencies, R_halo, M_dust, z):
        """
        (2darray, 1darray, [flex], [flex], [flex], [flex]) -> void

        [flex] = number / 1darray

        Updates the data in this instance and primes it for
        recalculation. Any data currently stored in this instance
        will be lost unless it has been copied elsewhere.

        L_nu: specific luminosity at all frequencies provided [ergs / s / Hz]
        frequencies: in [Hz]
        R_halo: radius of dark matter halo [cm]
        M_dust: dust mass of the galaxy(ies) [solar masses]
        z: redshift(s)

        All data and frequencies must be provided in ascending order of
        frequencies. The second dimension of the L_nu array can be any
        other property (such as different times, different galaxies, etc.)
        """

        DustModel.validate_data(L_nu, frequencies, R_halo, M_dust, z)

        self._L_nu = L_nu
        self._frequencies = frequencies
        self._R_halo = R_halo
        self._M_dust = M_dust
        self._z = z

        for key in self._modified.keys():
            self._modified[key] = True

    @property
    def kappa_nu(self):
        """
        (void) -> 2darray

        Calculates the dust opacity given the frequencies. The data
        returned is in cm^2 / g.
        """
        if self.must_recalculate('_kappa_nu'):
            self._kappa_nu = np.outer(0.1 * (self._frequencies / 1e9 / 1000)**2, np.ones(self._L_nu.shape[1]))
            self._modified['_kappa_nu'] = False
        return self._kappa_nu

    # TODO make this an input / get from parameter file
    @property
    def R_dust(self):
        """
        (void) -> num / 1darray

        Calculates and/or returns the dust radius / radii.
        """
        # TODO Check if ARES provides the galaxy radii already
        if self.must_recalculate('_R_dust'):
            self._R_dust = 0.018 * self.R_halo
            self._modified['_R_dust'] = False
        return self._R_dust

    @property
    def tau_nu(self):
        """
        (void) -> 2darray

        Calculates the optical depth of the dust given the
        stellar mass, metallicities and redshifts.
        """
        if self.must_recalculate('_tau_nu'):
            self._tau_nu = 3 * self.M_dust * g_per_msun / 4 / np.pi / self.R_dust**2
            self._tau_nu = np.outer(np.ones(self.L_nu.shape[0]), self._tau_nu) * self.kappa_nu
            self._modified['_tau_nu'] = False
        return self._tau_nu

    @property
    def f_geom(self):
        """
        (void) -> 2darray

        Calculates and/or returns the geometric factor associated with a
        homogeneously distributed galaxy.
        """
        if self.must_recalculate('_f_geom'):
            self._f_geom = (1 - np.exp(-self.tau_nu)) / self.tau_nu
            self._modified['_f_geom'] = False
        return self._f_geom

    @property
    def T_cmb(self):
        """
        (void) -> num / 1darray

        Calculates and returns the CMB temperature given redshift.
        """
        if self.must_recalculate('_T_cmb'):
            self._T_cmb = 2.725 * (1 + self.z)
            self._modified['_T_cmb'] = False
        return self._T_cmb

    @property
    def absorbed_power_per_dust_mass(self):
        """
        (void) -> 1darray

        Calculates and returns the total absorbed power by the dust.
        The units are ergs / s.
        """
        if self.must_recalculate('_abs_pow'):

            # non-cmb and cmb are temporary variables which are then added
            # and integrated to get the full power in ergs / s
            non_cmb = self.L_nu * self.f_geom * F_STAR * self.kappa_nu
            if isinstance(self.R_dust, np.ndarray):
                non_cmb /= np.outer(np.ones(self.L_nu.shape[0]), self.R_dust**2)
            else:
                non_cmb /= self.R_dust**2

            cmb = 8 * np.pi * h / c**2 * self.kappa_nu
            cmb *= np.outer(self.frequencies**3, np.ones(self.L_nu.shape[1]))
            cmb /= np.exp(h / k_B * np.outer(self.frequencies, 1 / self.T_cmb)) - 1

            self._abs_pow = scipy.integrate.simps(cmb + non_cmb, self.frequencies, axis = 0)
            self._modified['_abs_pow'] = False

        return self._abs_pow

    @property
    def T_dust(self):
        """
        (void) -> 1darray

        Calculates and returns the dust temperature from the total power.
        
        Does NOT unprime the instance for recalculation, the
        'emission_spectrum' property must be requested for that.
        """
        if self.must_recalculate('_T_dust'):
            prefactor = 64e-25 / 63 * np.pi**7 * k_B**6 / c**2 / h**5
            self._T_dust = (self.absorbed_power_per_dust_mass / prefactor)**(1/6)
            self._modified['_T_dust'] = False
        return self._T_dust

    def get_emission_band(self, nu_min, nu_max, num_samples):
        """
        (float, float, int) -> 1darray, 2darray

        Calculates and returns the emission spectrum in ergs / s / Hz
        in a given band with specified precision.

        The band must be provided in Hz.

        Does NOT unprime the instance for recalculation.
        """
        freqs = np.linspace(nu_min, nu_max, num_samples)
        kappa_nu_temp = 0.1 * (freqs / 1e9 / 1000)**2

        emissions = 8 * np.pi * h / c**2 * kappa_nu_temp * freqs**3
        emissions = np.outer(emissions, np.ones(self.L_nu.shape[1]))
        emissions /= np.exp(h / k_B * np.outer(freqs, 1 / self.T_dust)) - 1

        if isinstance(self.M_dust, np.ndarray):
                emissions *= np.outer(np.ones(self.L_nu.shape[0]), self.M_dust) * g_per_msun
        else:
            emissions *= self.M_dust * g_per_msun

        return freqs, emissions


    # TODO Fix this so we can get a good spectrum by default
    @property
    def emission_frequencies(self):
        """
        (void) -> 2darray

        Shifts the frequencies to the right so that the emission
        spectrum is fully within the range
        """
        if self.must_recalculate('_emission_freqs'):
            nu_peak = 5.879e10 * self.T_dust
            self._emission_freqs = np.linspace(1, nu_peak * 5, self.frequencies.shape[0])
            self._modified['_emission_freqs'] = False
        return self._emission_freqs

    @property
    def emission_spectrum(self):
        """
        (void) -> 1darray, 2darray
        
        Calculates and returns the emission spectrum of the dust in each galaxy
        in ergs / s / Hz.
        """
        if self.must_recalculate('_emission_spec'):

            kappa_nu_temp = 0.1 * (self.emission_frequencies / 1e9 / 1000)**2

            self._emission_spec = 8 * np.pi * h / c**2 * kappa_nu_temp * self.emission_frequencies**3
            self._emission_spec /= np.exp(h / k_B * self.emission_frequencies * np.outer(self.L_nu.shape[0], 1 / self.T_dust)) - 1

            if isinstance(self.M_dust, np.ndarray):
                self._emission_spec *= np.outer(np.ones(self.L_nu.shape[0]), self.M_dust) * g_per_msun
            else:
                self._emission_spec *= self.M_dust * g_per_msun

            self._modified['_emission_spec'] = False

        return self._emission_spec