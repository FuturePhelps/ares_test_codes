import ares
import numpy as np
import matplotlib.pyplot as plt
from ares.physics.Constants import c

pars = ares.util.ParameterBundle('mirocha2020:univ')
pop = ares.populations.GalaxyPopulation(**pars)

x1, phi1 = pop.LuminosityFunction(6, None, wave = 1600)
x2, phi2 = pop.LuminosityFunction(6, None, wave = 3e5)

plt.semilogy(x1, phi1, label = 'stellar luminosity (1600 A)')
plt.semilogy(x2, phi2, label = 'dust luminosity (300 000 A)')
plt.legend()
plt.show()