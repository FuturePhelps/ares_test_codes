import ares
from DustPopulation import DustPopulation
import matplotlib.pyplot as plt
import numpy as np

pars = ares.util.ParameterBundle('mirocha2020:univ')

# For other parameters
pop = ares.populations.GalaxyPopulation(**pars)

dust = DustPopulation(pop)

freqs, SED = dust.dust_sed(1e10, 3e13, 1000)

fig, ax = plt.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Spectral Radiance (ergs / s / Hz)")
ax.set_title("SED over redshift for one galaxy\nin the ensemble")

for i in range(len(dust.z)):
    ax.loglog(freqs, SED[10,:,i], label = ('z = %.1f' % dust.z[i]))
ax.legend()
fig.savefig("figures/SED_some_zero_dust_mass_log.png", bbox_inches = 'tight')
fig.savefig("figures/SED_some_zero_dust_mass_log.pdf", bbox_inches = 'tight')