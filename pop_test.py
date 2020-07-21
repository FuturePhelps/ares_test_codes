import ares
import numpy as np
import matplotlib.pyplot as pl
from ares.physics.Constants import c

pars = ares.util.ParameterBundle('mirocha2020:univ')

# For other parameters
pop = ares.populations.GalaxyPopulation(**pars)

zobs = 9.5
freqs = np.geomspace(1e14, 1e17, 900)
# This has to be a wavelength array
waves = (c / freqs) * 1e8

# Just need SFRs and times (or redshifts)
hist = pop.histories

# Compute the spectra for all objects
# This is already in ergs / s / Hz
spectra = pop.synth.Spectrum(waves, zobs=zobs, sfh=hist['SFR'], tarr=hist['t'])

print(spectra.shape)

fig, ax = pl.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Spectral Radiance (ergs / s / Hz)")
ax.set_title("Galactic spectrum from ares.static.SpectralSynthesis")
ax.loglog(freqs, spectra[-1])
fig.savefig('figures/SpectralSynthesis.pdf', bbox_inches = 'tight')
fig.savefig('figures/SpectralSynthesis.png', bbox_inches = 'tight')

Ms = pop.get_field(zobs, 'Ms')

print(Ms[-1] / 1e10, Ms[-100] / 1e10)