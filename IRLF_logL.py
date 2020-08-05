import ares
import numpy as np
import matplotlib.pyplot as plt
from ares.physics.Constants import c

pars = ares.util.ParameterBundle('mirocha2020:univ')
pop = ares.populations.GalaxyPopulation(**pars)

x, phi = pop.LuminosityFunction(4, None, wave = 2.5e6, mags = False)

fig, ax = plt.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.semilogy(x, phi, label = r'dust luminosity (250 $\mu\mathrm{m}$)')
ax.legend()
ax.set_title(r"Galaxy Luminosity Functions")
ax.set_xlabel(r"$\log{\left(L/L_\odot\right)}$")
ax.set_ylabel(r"$\phi$ [$\mathrm{Mpc}^{-3}\mathrm{dex}^{-1}$]")
fig.savefig("../figures/IRLF_LogL.pdf", bbox_inches = 'tight')
fig.savefig("../figures/IRLF_LogL.png", bbox_inches = 'tight')