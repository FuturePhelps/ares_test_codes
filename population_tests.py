import ares
import numpy as np
import matplotlib.pyplot as plt
import os
import DustModel
from ares.physics.Constants import g_per_msun, cm_per_kpc

pars = ares.util.ParameterBundle('mirocha2017:base').pars_by_pop(0, 1)

pars['pop_sed'] = 'eldridge2009'
pars['pop_Z'] = 0.02

pop = ares.populations.GalaxyPopulation(**pars)
src = pop.src # will be the same as in previous example

data = (src.data[:,0] * src.dwdn)[::-1]
freqs = np.sort(src.frequencies)

plt.plot(freqs, data)
plt.xlim(0,2e16)
plt.show()