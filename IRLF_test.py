import ares
import numpy as np
import matplotlib.pyplot as plt

pars = ares.util.ParameterBundle('mirocha2020:univ')
pop = ares.populations.GalaxyPopulation(**pars)

dust = pop.dust

dust.L_nu

for_fun_1 = dust.Luminosity(4, 1500)