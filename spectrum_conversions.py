import ares
import numpy as np
import matplotlib.pyplot as plt

src = ares.sources.SynthesisModel(source_sed='eldridge2009', source_Z=0.02)

# src.wavelengths: increasing order, angstroms
# src.frequencies: decreasing order, Hz

power_per_wavelength_per_SFR = src.data # power / angstrom / SFR
power_per_frequency_per_SFR = power_per_wavelength_per_SFR * np.outer(src.dwdn, np.ones(src.data.shape[1])) # power / Hz / SFR

plt.loglog(src.wavelengths, power_per_wavelength_per_SFR[:,0], label = 'wavelength')
plt.legend()
plt.show()

plt.loglog(src.frequencies[::-1], power_per_frequency_per_SFR[::-1,0], label = 'frequency')
plt.title(r"Specific Luminosity $L_\nu$ vs frequency")
plt.xlabel("Frequency (Hz)")
plt.ylabel(r"$L_\nu$ (ergs / s / Hz)")
plt.legend()
plt.show()