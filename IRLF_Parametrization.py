import ares
import numpy as np
import matplotlib.pyplot as plt

z = np.arange(4, 9, 1)

pars = ares.util.ParameterBundle('mirocha2020:univ')
pars['pop_dust_zmin'] = z[0]
pars['pop_dust_zmax'] = z[-1]
pars['pop_dust_Nz'] = len(z)

pop1 = ares.populations.GalaxyPopulation(**pars)
dust1 = pop1.dust

pars['pop_dust_experimental'] = True

pop2 = ares.populations.GalaxyPopulation(**pars)
dust2 = pop2.dust

x1, phi1 = pop1.LuminosityFunction(4, None, wave = 2.5e6, mags = False)
x2, phi2 = pop2.LuminosityFunction(4, None, wave = 2.5e6, mags = False)

fig, ax = plt.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.plot(x1, np.log10(phi1), label = r'dust luminosity ($z=4$), Imara et al. (2018)')
ax.plot(x2, np.log10(phi2), label = r'dust luminosity ($z=4$), Log-linear Parametrization')
ax.legend()
ax.set_title(r"Dust Luminosity Functions")
ax.set_xlabel(r"$\log{\left(L_{250}/L_\odot\right)}$")
ax.set_ylabel(r"$\log{\left(\phi\right)}$ [$\mathrm{Mpc}^{-3}\mathrm{dex}^{-1}$]")
fig.savefig("../figures/IRLF_Parametrization.pdf", bbox_inches = 'tight')
fig.savefig("../figures/IRLF_Parametrization.png", bbox_inches = 'tight')