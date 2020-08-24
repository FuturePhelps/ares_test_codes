import ares
import numpy as np
import matplotlib.pyplot as plt
from ares.physics.Constants import h, k_B, c, Lsun

z = np.arange(4, 9, 1)

pars = ares.util.ParameterBundle('mirocha2020:univ')
pars['pop_dust_zmin'] = z[0]
pars['pop_dust_zmax'] = z[-1]
pars['pop_dust_Nz'] = len(z)
pars['pop_dust_experimental'] = True
pop = ares.populations.GalaxyPopulation(**pars)
dust = pop.dust


M_star = []
for i in z:
    M_star.append(pop.get_field(i, 'Ms'))
M_star = np.array(M_star).transpose()

T_dust = dust.T_dust

fig, ax = plt.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
for i in range(len(z)):
    ax.semilogx(M_star[:,i], T_dust[:,i], '.', label = r'$z=\;$' + ('%d' % z[i]))
ax.legend()
ax.set_title(r"Dust Temperature vs Stellar Mass")
ax.set_xlabel(r"$M_{\mathrm{star}}\;\;\left[M_\odot\right]$")
ax.set_ylabel(r"$T_{\mathrm{dust}}\;\;\left[\mathrm{K}\right]$")
fig.savefig("../figures/Tdust_vs_Mstar.pdf", bbox_inches = 'tight')
fig.savefig("../figures/Tdust_vs_Mstar.png", bbox_inches = 'tight')
plt.clf()


SFR = np.zeros((dust.Ngalaxies, dust.Nz))
for i in range(dust.Nz):
    SFR[:,i] = pop.get_field(z[i], 'SFR')

np.savetxt('../test_files/DGR.csv', dust.DGR, '%.8f', '\t')
np.savetxt('../test_files/Z.csv', dust.Z, '%.8f', '\t')
np.savetxt('../test_files/SFR.csv', SFR, '%.8f', '\t')
np.savetxt('../test_files/tau_nu.csv', dust.tau_nu[-1,:,:], '%.8f', '\t')
np.savetxt('../test_files/f_geom.csv', (1 - np.exp(-dust.tau_nu[-1,:,:])) / dust.tau_nu[-1,:,:], '%.8f', '\t')

print("M_star =", M_star[15:,:])
print("Dust mass =", dust.M_dust[15:,:])
print("DGR =", dust.DGR)
print("Z =", dust.Z)

fig, ax = plt.subplots(1,1)
fig.set_size_inches(6.5, 3.7)
for i in range(len(z)):
    ax.plot(M_star[:,i], dust.M_dust[:,i], '.', label = r'$z=\;$' + ('%d' % z[i]))
ax.legend()
ax.set_title(r"Dust Mass vs Stellar Mass")
ax.set_xlabel(r"$M_{\mathrm{star}}\;\;\left[M_\odot\right]$")
ax.set_ylabel(r"$M_{\mathrm{dust}}\;\;\left[M_\odot\right]$")
fig.savefig("../figures/Mdust_vs_Mstar.pdf", bbox_inches = 'tight')
fig.savefig("../figures/Mdust_vs_Mstar.png", bbox_inches = 'tight')
plt.clf()