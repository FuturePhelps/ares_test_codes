# This file always changes, it is only used to probe ARES
# properties and available data

import ares
import numpy as np
import matplotlib.pyplot as plt

pars = ares.util.ParameterBundle('mirocha2020:univ')
pop = ares.populations.GalaxyPopulation(**pars)\

M_star = np.zeros((pop.dust.Ngalaxies, pop.dust.Nz))
for i in range(pop.dust.Nz):
    M_star[:,i] = pop.get_field(pop.dust.z[i], 'Ms')

plt.semilogx(M_star[:,-1], pop.dust.T_dust[:,-1], '.')
plt.show()