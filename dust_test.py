import ares
import matplotlib.pyplot as pl

pars = ares.util.ParameterBundle('mirocha2020:univ')

pop = ares.populations.GalaxyPopulation(**pars)

dust = pop.dust

freqs, seds = pop.dust_sed(1e10, 3e13, 1000)

for i in range(len(dust.z)):
    pl.loglog(freqs, seds[-1,:,i], label = ("z = %.1f" % dust.z[i]))
pl.legend()
pl.show()