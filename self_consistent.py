import ares
import numpy as np
import matplotlib.pyplot as pl

pars = ares.util.ParameterBundle('mirocha2020:univ')

pop = ares.populations.GalaxyPopulation(**pars)

# First, some data
obslf = ares.analysis.GalaxyPopulation()
ax = obslf.PlotLF(z=6, round_z=0.2)

# Now, the predicted/calibrated UVLF
mags = np.arange(-25, -5, 0.1)
phi = pop.LuminosityFunction(6, mags)

ax.semilogy(mags, phi)

axes = obslf.PlotColors(pop, fig=2)

mags = pop.Magnitude(6., presets='hst')
beta = pop.Beta(6., presets='hst')

fig3, ax3 = pl.subplots(1, 1, num=3)
ax3.scatter(mags, beta, alpha=0.1, color='b', edgecolors='none')

mags = np.arange(-25, -10, 0.5) # bin centers
beta, beta_s = pop.Beta(6., presets='hst', Mbins=mags, return_binned=True,
        return_scatter=True)

# Plot scatter in each MUV bin as errorbars
ax3.errorbar(mags, beta, yerr=beta_s.T, color='b', marker='s', fmt='o')
pl.show()