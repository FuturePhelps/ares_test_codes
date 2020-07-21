from DustModel import DustModel as DM
import ares
from ares.physics.Constants import g_per_msun, cm_per_kpc
import numpy as np
import matplotlib.pyplot as plt
import os

path = os.path.realpath(os.path.dirname(__file__))
os.chdir(path)

pars = ares.util.ParameterBundle('mirocha2017:base').pars_by_pop(0, 1)

pars['pop_sed'] = 'eldridge2009'
pars['pop_Z'] = 0.007

pop = ares.populations.GalaxyPopulation(**pars)
src = pop.src # will be the same as in previous example

SFR = 160

L_nu = (src.data * np.outer(src.dwdn, np.ones(src.data.shape[1])) * SFR)[::-1]
frequencies = src.frequencies[::-1]
R_halo = 0.57 * cm_per_kpc / 0.018
M_star = 7.9e9
z = 9.5
Z = 0.007

dust = DM(L_nu, frequencies, R_halo, M_star, z, Z)

R_gal = dust.R_dust / cm_per_kpc
print("M_star =", dust.M_star / 1e10)
print("SFR =", SFR)
print("R_gal =", R_gal)
print("Z =", dust.Z)
print("DGR =", dust.DGR * 463)
print("M_gas =", dust.M_gas / 1e10)
print("M_dust =", dust.M_dust / 1e7)
print("T_dust =", dust.T_dust[0])

freqs, emissions = dust.get_emission_band(1e10, 1e15, 50000)

fig, ax = plt.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.set_title("Dust emissions and stellar emissions")
ax.loglog(freqs, emissions[:,0], 'r-', label = "Dust Emissions")
ax.loglog(frequencies, L_nu[:,0], 'b-', label = "Stellar Emissions")
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Spectral Radiance (ergs / s / Hz)")
ax.legend()
fig.savefig("figures/dust_emissions_2020-07-03.pdf", bbox_inches = 'tight')
fig.savefig("figures/dust_emissions_2020-07-03.png", bbox_inches = 'tight')

fig, ax = plt.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.set_title("Dust emissions and stellar emissions")
ax.set_ylim(np.max(emissions[:,0]) / 1000, np.max(emissions[:,0]) * 100)
ax.loglog(freqs, emissions[:,0], 'r-', label = "Dust Emissions")
ax.loglog(frequencies, L_nu[:,0], 'b-', label = "Stellar Emissions")
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Spectral Radiance (ergs / s / Hz)")
ax.legend()
fig.savefig("figures/dust_emissions_tight_2020-07-03.pdf", bbox_inches = 'tight')
fig.savefig("figures/dust_emissions_tight_2020-07-03.png", bbox_inches = 'tight')