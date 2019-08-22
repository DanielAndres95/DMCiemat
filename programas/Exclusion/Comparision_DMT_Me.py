import numpy as np
import matplotlib.pyplot as plt
import math as math
import global_functions as glb

sec = np.asarray(glb.read_datatxt("secdmtools.txt", " "))
mx = np.asarray(glb.read_datatxt("massdmtools.txt", " "))
print(mx)

Z = 18 
N = 22
A = Z + N
M = 931.5*A*1e-3 #Masa del nucleo en GeVs
GF = 1.166e-11*1e6 #Constante de Fermi en GeVs
mr = (mx*M/(mx + M)) #Masa reducida
mrn = (mx*(M/A) /(mx + M/A))
r = 4*mr**2/(mx*M)
rho = 0.3 # GeV cm⁻³
Ro = 3000 # event / kg / day
vo = 230 * 1e5 *3600 * 24 #cm day⁻¹
No = 6.022e23*1e3 / A # kg⁻¹

# Total cross section #
sec_Napprox = (math.sqrt(math.pi) * mx * Ro / (2*No*rho*vo*A)) * (mrn / mr)**2

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))	

ax.plot(mx, sec, 'lime',lw=3)
ax.plot(mx, sec_Napprox*1e28)
ax.grid('--', alpha=0.3)

plt.ylabel(r"$\sigma\ (cm^2)$",fontsize=14)
plt.xlabel(r'$m_\chi\ (GeV)$',fontsize=14)
plt.title(r'$\rho = 0.3\ (GeV\cdot cm^{-3}),\ R_0 = 3000\ (events\cdot day^{-1}\cdot kg^{-1}) $')

plt.legend([r'DM Tools', r'Theoretical'])
plt.xscale('log')
plt.yscale('log')

plt.show()