import ares
import numpy as np
import matplotlib.pyplot as plt
from ares.physics.Constants import c

pars = ares.util.ParameterBundle('mirocha2020:ereg_eduty_edtmr')
pars['pop_dust_experimental'] = True
#pars['pop_dust_distrib'] = 'pt src'
#pars['pop_dust_yield'] = 40
#pars['pq_func_par3[0]'] = -0.01
pop = ares.populations.GalaxyPopulation(**pars)

x, phi = pop.LuminosityFunction(4, None, wave = 2.5e6, mags = False)

data_x = np.array([10, 10.25, 10.50, 10.75])
data_y = np.array([-3.52, -3.36, -3.68, -4.10])
errors = np.array([[0.76, 0.47, 0.45, 0.76],[0.59, 0.48, 0.48, 0.68]])

fig, ax = plt.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.plot(x, np.log10(phi), label = r'dust luminosity ($z=4$)')
ax.errorbar(data_x, data_y, errors, None, 'ro', 'k', capsize = 2, label = r"ALPINE data ($3.5\leq z\leq 4.5$)")
ax.legend()
ax.set_title(r"Dust Luminosity Function and ALPINE data")
ax.set_xlabel(r"$\log{\left(L_{250}/L_\odot\right)}$")
ax.set_ylabel(r"$\log{\left(\phi\right)}$ [$\mathrm{Mpc}^{-3}\mathrm{dex}^{-1}$]")
fig.savefig("../figures/IRLF_ALPINE_z4.pdf", bbox_inches = 'tight')
fig.savefig("../figures/IRLF_ALPINE_z4.png", bbox_inches = 'tight')


x, phi = pop.LuminosityFunction(5.25, None, wave = 2.5e6, mags = False)

data_x = np.array([10.25, 10.50, 10.75, 11.00, 11.25, 11.50, 11.75])
data_y = np.array([-3.49, -3.21, -3.61, -3.52, -3.61, -3.91, -3.91])
errors = np.array([[0.65, 0.37, 0.44, 0.58, 0.47, 0.62, 0.62],[0.45, 0.30, 0.47, 0.41, 0.34, 0.46, 0.46]])

fig, ax = plt.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.plot(x, np.log10(phi), label = r'dust luminosity ($z=5.25$)')
ax.errorbar(data_x, data_y, errors, None, 'ro', 'k', capsize = 2, label = r"ALPINE data ($4.5\leq z\leq 6.0$) [All emitters]")
ax.legend()
ax.set_title(r"Dust Luminosity Function and ALPINE data")
ax.set_xlabel(r"$\log{\left(L_{250}/L_\odot\right)}$")
ax.set_ylabel(r"$\log{\left(\phi\right)}$ [$\mathrm{Mpc}^{-3}\mathrm{dex}^{-1}$]")
fig.savefig("../figures/IRLF_ALPINE_z525_all.pdf", bbox_inches = 'tight')
fig.savefig("../figures/IRLF_ALPINE_z525_all.png", bbox_inches = 'tight')

x, phi = pop.LuminosityFunction(5.25, None, wave = 2.5e6, mags = False)

data_x = np.array([10.25, 10.50, 10.75, 11.00, 11.25])
data_y = np.array([-3.45, -3.31, -3.62, -3.51, -3.70])
errors = np.array([[0.58, 0.51, 0.44, 0.59, 0.45],[0.47, 0.41, 0.47, 0.41, 0.48]])

fig, ax = plt.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.plot(x, np.log10(phi), label = r'dust luminosity ($z=5.25$)')
ax.errorbar(data_x, data_y, errors, None, 'ro', 'k', capsize = 2, label = r"ALPINE data ($4.5\leq z\leq 6.0$) [C II emitters only]")
ax.legend()
ax.set_title(r"Dust Luminosity Function and ALPINE data")
ax.set_xlabel(r"$\log{\left(L_{250}/L_\odot\right)}$")
ax.set_ylabel(r"$\log{\left(\phi\right)}$ [$\mathrm{Mpc}^{-3}\mathrm{dex}^{-1}$]")
fig.savefig("../figures/IRLF_ALPINE_z525_cii.pdf", bbox_inches = 'tight')
fig.savefig("../figures/IRLF_ALPINE_z525_cii.png", bbox_inches = 'tight')

x, phi = pop.LuminosityFunction(4, None, wave = 2.5e6, mags = False, total_IR = True)

data_x = np.array([11.75, 12.00, 12.25, 12.50])
data_y = np.array([-3.37, -3.37, -4.10, -4.10])
errors = np.array([[0.82, 0.58, 0.78, 0.78],[0.40, 0.40, 0.59, 0.52]])

fig, ax = plt.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.plot(x, np.log10(phi), label = r'total dust luminosity ($z=4$)')
ax.errorbar(data_x, data_y, errors, None, 'ro', 'k', capsize = 2, label = r"ALPINE data ($3.5\leq z\leq 4.5$)")
ax.legend()
ax.set_title(r"Total Dust Luminosity Function and ALPINE data")
ax.set_xlabel(r"$\log{\left(L_{IR}/L_\odot\right)}$")
ax.set_ylabel(r"$\log{\left(\phi\right)}$ [$\mathrm{Mpc}^{-3}\mathrm{dex}^{-1}$]")
fig.savefig("../figures/IRLF_ALPINE_total_z4.pdf", bbox_inches = 'tight')
fig.savefig("../figures/IRLF_ALPINE_total_z4.png", bbox_inches = 'tight')

x, phi = pop.LuminosityFunction(5.25, None, wave = 2.5e6, mags = False, total_IR = True)

data_x = np.array([11.25, 11.50, 11.75, 12.00, 12.25, 12.50])
data_y = np.array([-3.91, -3.91, -3.45, -3.34, -3.75, -4.12])
errors = np.array([[0.78, 0.78, 0.53, 0.40, 0.52, 0.78],[0.55, 0.54, 0.39, 0.32, 0.38, 0.67]])

fig, ax = plt.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.plot(x, np.log10(phi), label = r'total dust luminosity ($z=5.25$)')
ax.errorbar(data_x, data_y, errors, None, 'ro', 'k', capsize = 2, label = r"ALPINE data ($4.5\leq z\leq 6.0$) [C II emitters only]")
ax.legend()
ax.set_title(r"Total Dust Luminosity Function and ALPINE data")
ax.set_xlabel(r"$\log{\left(L_{IR}/L_\odot\right)}$")
ax.set_ylabel(r"$\log{\left(\phi\right)}$ [$\mathrm{Mpc}^{-3}\mathrm{dex}^{-1}$]")
fig.savefig("../figures/IRLF_ALPINE_total_z525_cii.pdf", bbox_inches = 'tight')
fig.savefig("../figures/IRLF_ALPINE_total_z525_cii.png", bbox_inches = 'tight')

x, phi = pop.LuminosityFunction(5.25, None, wave = 2.5e6, mags = False, total_IR = True)

data_x = np.array([11.25, 11.50, 11.75, 12.00, 12.25, 12.50, 12.75, 13.00, 13.25])
data_y = np.array([-3.94, -3.91, -3.45, -3.38, -3.79, -4.16, -3.76, -3.63, -4.21])
errors = np.array([[0.78, 0.78, 0.53, 0.48, 0.63, 0.78, 0.78, 0.79, 0.78],[0.55, 0.54, 0.38, 0.34, 0.41, 0.59, 0.55, 0.43, 0.55]])

fig, ax = plt.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.plot(x, np.log10(phi), label = r'total dust luminosity ($z=5.25$)')
ax.errorbar(data_x, data_y, errors, None, 'ro', 'k', capsize = 2, label = r"ALPINE data ($4.5\leq z\leq 6.0$) [All emitters]")
ax.legend()
ax.set_title(r"Total Dust Luminosity Function and ALPINE data")
ax.set_xlabel(r"$\log{\left(L_{IR}/L_\odot\right)}$")
ax.set_ylabel(r"$\log{\left(\phi\right)}$ [$\mathrm{Mpc}^{-3}\mathrm{dex}^{-1}$]")
fig.savefig("../figures/IRLF_ALPINE_total_z525_all.pdf", bbox_inches = 'tight')
fig.savefig("../figures/IRLF_ALPINE_total_z525_all.png", bbox_inches = 'tight')