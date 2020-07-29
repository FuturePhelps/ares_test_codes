import ares
import numpy as np
import matplotlib.pyplot as plt
from ares.physics.Constants import c

pars = ares.util.ParameterBundle('mirocha2020:univ')
pop = ares.populations.GalaxyPopulation(**pars)

freqs, sed = pop.dust_sed(1e10, 1e15, 1000)

plt.loglog(freqs, sed[-1,:,2])
plt.loglog(pop.dust.frequencies, pop.dust.L_nu[-1,:,2])
plt.show()

wave_max = pop.src.wavelengths.max()
wave = c / 1e13 * 1e8
L = pop.Luminosity(6, wave)
pop.dust.L_nu
L_dust = pop.dust.Luminosity(6, wave)
print(wave_max)
print(L)
print(L_dust)
print(L_dust > L)