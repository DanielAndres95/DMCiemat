import numpy as np
import matplotlib.pyplot as plt
#from basic_units import radians
import math as mth

# ER(angle):

theta = np.linspace(0, mth.pi, 100)
vo = 0.001 # v/c
M = 39.96 * 931.5 # 40Ar mass in Mev's/c2
m = 100 * 1e3 # 100 Gev's/c2 WIMP mass in Mev's/c2

mr = m*M/(m + M) # Reduced mass
ER = (mr**2 * vo**2 * (1 - np.cos(theta)) ) / M # Non-relativistic recoil energy

plt.figure(figsize=(8, 6))
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.tight_layout()
plt.plot(theta, ER, 'chocolate')
plt.grid(linestyle='--', alpha=0.4)
plt.xlabel(r'$\theta$', fontsize=14)# xunits=radians)
plt.ylabel(r'$E_R (MeV)$', fontsize=14)
plt.text(2.6, 0.002, r'$m_\chi = 100 GeV$', fontsize=14)
plt.savefig('/afs/ciemat.es/user/t/tomas/CIEMAT/Dark_Matter/programas/figures/ER(theta).pdf')

# ER(m):
m = np.linspace(10, 1000, 100) * 1e3 # 10 GeV - 1 TeV

mr = m*M/(m + M) # Reduced mass
theta  = mth.pi
ER = (mr**2 * vo**2 * (1 - np.cos(theta)) ) / M # Non-relativistic recoil energy

plt.figure(figsize=(8, 6))
plt.plot(m/1e3, ER*1e3, 'darkblue')
plt.grid(linestyle='--', alpha=0.4)
plt.xlabel(r'$m_\chi (GeV)$', fontsize=14)
plt.ylabel(r'$E_R (keV)$', fontsize=14)
plt.text(9e5, 0.005, r'$\theta = \pi$',
         bbox=dict(facecolor='white', edgecolor='white', alpha=1), fontsize=14)
plt.savefig('/afs/ciemat.es/user/t/tomas/CIEMAT/Dark_Matter/programas/figures/ER(mx).pdf')

plt.show()


 

