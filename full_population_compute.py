import ares
from DustPopulation import DustPopulation
import matplotlib.pyplot as plt
import numpy as np

# Create a GalaxyEnsemble object from the parameters
pars = ares.util.ParameterBundle('mirocha2020:univ')
pop = ares.populations.GalaxyPopulation(**pars)

# This is how we access the DustPopulation object
dust = pop.dust

# Retrieve the dust SEDs between nu = 1e10 and 3e13 with 1000 sampled frequencies
freqs, SED = pop.dust_sed(1e10, 3e13, 1000)

# Plot the seds for one galaxy at all sampled redshifts
fig, ax = plt.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Spectral Radiance (ergs / s / Hz)")
ax.set_title("SED over redshift for one galaxy\nin the ensemble")

for i in range(dust.Nz):
    ax.loglog(freqs, SED[-1,:,i], label = ('z = %.1f' % dust.z[i]))
ax.legend()

# Save it in the figures folder as both pdf and png
fig.savefig("figures/SED_galactic_dust.png", bbox_inches = 'tight')
fig.savefig("figures/SED_galactic_dust.pdf", bbox_inches = 'tight')