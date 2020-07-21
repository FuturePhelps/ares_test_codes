import ares
import numpy as np
import matplotlib.pyplot as pl
import os
import scipy
from ares.physics.Constants import c, h, k_B, g_per_msun

path = os.path.realpath(os.path.dirname(__file__))
os.chdir(path)

src = ares.sources.SynthesisModel(source_sed='eldridge2009', source_Z=0.02)

pl.figure(1)
pl.title("Galaxy SED vs wavelength at different times")
pl.xlabel(r"Wavelength ($\AA$)")
pl.ylabel(r"SED (ergs / s / $\AA$ / ($M_\odot$ / yr)")
for i in range(3):
    t = src.times[i]
    pl.loglog(src.wavelengths, src.data[:,i], label=r'$t = {}$ Myr'.format(t))

pl.legend()
pl.savefig("figures/example_1.pdf", bbox_inches = 'tight')
pl.savefig("figures/example_1.png", bbox_inches = 'tight')

pl.figure(2)
pl.title("Galaxy SED vs time at different wavelengths")
pl.xlabel(r"Time (Myr)")
pl.ylabel(r"SED (ergs / s / $\AA$ / ($M_\odot$ / yr)")
for i in range(0, 1000, 200):
    wave = src.wavelengths[i]
    pl.loglog(src.times, src.data[i,:], label=r'$\lambda = {} \AA$'.format(wave))

pl.legend()
pl.savefig("figures/example_2.pdf", bbox_inches = 'tight')
pl.savefig("figures/example_2.png", bbox_inches = 'tight')

data = src.data[:,0] # Let's just look at one time for now

frequencies = src.frequencies[::-1]

# Estimating dust radius
R_halo = 3.086e+23 # cm
R_dust = 0.018 * R_halo # cm

# Picking whatever value for redshift
z = 9.5

# pick some constant SFR for now
SFR = 1 # solar mass per year

M_star = 7.9e9 # Solar masses
M_gas = 3.87e9 * (1 + z)**1.35 * (M_star / 1e10)**0.49 # solar masses
DGR = 1 / 463
M_dust = M_gas * DGR # solar masses

kappa_nu = 0.1 * (frequencies / 1e9 / 1000)**2

fig, ax = pl.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.plot(frequencies, kappa_nu)
ax.set_title(r"Dust Opacity vs Frequency")
ax.set_xlabel(r"Frequency $\nu$ (Hz)")
ax.set_ylabel(r"Dust Opacity $\kappa_\nu$ ($\mathrm{cm}^2/\mathrm{g}$)")
fig.savefig("figures/kappa_nu_vs_nu.pdf", bbox_inches = 'tight')
fig.savefig("figures/kappa_nu_vs_nu.png", bbox_inches = 'tight')

tau_nu = 3 * (M_dust * g_per_msun) * kappa_nu / (4 * np.pi * R_dust**2)

fig, ax = pl.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.plot(frequencies, tau_nu)
ax.set_title(r"Optical Depth vs Frequency")
ax.set_xlabel(r"Frequency $\nu$ (Hz)")
ax.set_ylabel(r"Optical Depth $\tau_\nu$ (dimensionless)")
fig.savefig("figures/tau_nu_vs_nu.pdf", bbox_inches = 'tight')
fig.savefig("figures/tau_nu_vs_nu.png", bbox_inches = 'tight')

# Geometric factor for homogeneous galaxy
f_geom = (1 - np.exp(-tau_nu)) / tau_nu

fig, ax = pl.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.plot(frequencies, f_geom)
ax.set_title(r"Geometric Factor vs Frequency")
ax.set_xlabel(r"Frequency $\nu$ (Hz)")
ax.set_ylabel(r"Geometric Factor $f_\mathrm{geom}$ (dimensionless)")
fig.savefig("figures/f_geom_vs_nu.pdf", bbox_inches = 'tight')
fig.savefig("figures/f_geom_vs_nu.png", bbox_inches = 'tight')

# Covering factor of interstellar dust
f_star = 0.25

# steps:
# 1. convert data from power / wavelength to power / frequency using src.dwdn
# 2. multiply by constant SFR
# 3. get the data at the given nu

L_nu = (data * src.dwdn * SFR)[::-1]

fig, ax = pl.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.loglog(frequencies, L_nu)
ax.set_title(r"Specific Luminosity $L_\nu$ vs frequency")
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel(r"$L_\nu$ (ergs / s / Hz)")
fig.savefig("figures/L_nu_vs_f.pdf", bbox_inches = 'tight')
fig.savefig("figures/L_nu_vs_f.png", bbox_inches = 'tight')

# Now I guess we can calculate the function we need to integrate
non_cmb = L_nu * f_geom * f_star * kappa_nu / R_dust**2

fig, ax = pl.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.loglog(frequencies, non_cmb)
ax.set_title("Non-CMB Absorbtion vs Frequency")
ax.set_xlabel(r"Frequency $\nu$ (Hz)")
ax.set_ylabel("Non-CMB Absorbtion per frequency (ergs / s / g / Hz)")
fig.savefig("figures/non_cmb_absorb_vs_nu.pdf", bbox_inches = 'tight')
fig.savefig("figures/non_cmb_absorb_vs_nu.png", bbox_inches = 'tight')

# Let's now do the cmb part
T_cmb = 2.725 * (1 + z)
cmb = 8 * np.pi * h / c**2 * frequencies**3 / (np.exp(h * frequencies / (k_B * T_cmb)) - 1) * kappa_nu

fig, ax = pl.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.loglog(frequencies, cmb)
ax.set_title("CMB Absorbtion vs Frequency")
ax.set_xlabel(r"Frequency $\nu$ (Hz)")
ax.set_ylabel("CMB Absorbtion per frequency (ergs / s / g / Hz)")
fig.savefig("figures/cmb_absorb_vs_nu.pdf", bbox_inches = 'tight')
fig.savefig("figures/cmb_absorb_vs_nu.png", bbox_inches = 'tight')

# Add both
absorbtion = non_cmb + cmb
fig, ax = pl.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.loglog(frequencies, absorbtion)
ax.set_title("Absorbtion vs Frequency")
ax.set_xlabel(r"Frequency $\nu$ (Hz)")
ax.set_ylabel("Absorbtion per frequency (ergs / s / g / Hz)")
fig.savefig("figures/absorb_vs_nu.pdf", bbox_inches = 'tight')
fig.savefig("figures/absorb_vs_nu.png", bbox_inches = 'tight')

# Estimate total power per unit mass of dust
P_total_per_dust_mass = scipy.integrate.simps(absorbtion, frequencies) # ergs / s / g

# Now how do I extract the dust temp... okay I think I got it
# I can do the integral in Wolfram alpha because only kappa_nu
# depends on frequency like nu^2, the rest is all constants.
# So we get that P_tot_per_dust_mass = 1.8513e-7 * T_dust**6 ergs / s / g:

prefactor = 64e-25 / 63 * np.pi**7 * k_B**6 / c**2 / h**5
T_dust = (P_total_per_dust_mass / prefactor)**(1/6)
print(T_dust)

# Calculate the emission spectrum
emission = 8 * np.pi * h / c**2 * frequencies**3 / (np.exp(h * frequencies / k_B / T_dust) - 1) * kappa_nu
fig, ax = pl.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.loglog(frequencies, emission)
ax.set_title("Dust light emissions per dust mass")
ax.set_xlabel(r"Frequency $\nu$ (Hz)")
ax.set_ylabel(r"Emission spectrum (ergs / s / Hz / g)")
fig.savefig("figures/emission_spectrum.pdf", bbox_inches = 'tight')
fig.savefig("figures/emission_spectrum.png", bbox_inches = 'tight')

# Let's change the units to the same as the original SED and compare them
emission_lambda = emission[::-1] * M_dust * g_per_msun / src.dwdn / SFR

fig, ax = pl.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
ax.loglog(src.wavelengths, emission_lambda, 'r-', label = "Dust Emissions")
ax.loglog(src.wavelengths, data, 'b-', label = "ARES SED")
ax.set_title("Galaxy SED")
ax.set_xlabel(r"Wavelength ($\AA$)")
ax.set_ylabel(r"Emission spectrum (ergs / s / $\AA$ / ($M_\odot$ / yr))")
ax.legend()
fig.savefig("figures/total_seds.pdf", bbox_inches = 'tight')
fig.savefig("figures/total_seds.png", bbox_inches = 'tight')