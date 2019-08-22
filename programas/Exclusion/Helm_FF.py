import numpy as np
import matplotlib.pyplot as plt 
import math as math
import scipy.special as sp 

# Constants:
m = 1e3 # 1 TeV WIMP
A = [132, 32, 40, 23]
plt.figure(figsize=(8, 6))	
hc = 197.3e-3 # GeV fm
v = 0.005 # units of c
theta = np.linspace(0, math.pi, num=1e3)

colors = ['forestgreen', 'b', 'r', 'k']
i = 0

for A in A:

	M = 0.932 * A # Nuclear mass in GeVs

	# Kinematics
	mr = m*M/(m+M)
	E_R = mr**2 * v**2 * (1 - np.cos(theta)) / M 
	q = np.sqrt(2 * M * E_R) #GeV

	# Form Factor definitions:
	c = 1.23*A**(1/3) - 0.60 #fm
	s = 0.9 #fm
	a = 0.52 #fm
	rn = np.sqrt(c**2 + 7/3*math.pi**2*a**2 - 5*s**2)

	# Form factor:
	j1 = (np.sin(q*rn/hc) - q*rn/hc*np.cos(q*rn/hc)) / (q*rn/hc)**2
	F = 3 * j1 * np.exp(-(q*s/hc)**2/2) / (q*rn/hc)

	# Plot:
	color = colors[i]
	i += 1
	plt.rc('text',usetex=True)
	plt.rc('font',family='serif')
	plt.tight_layout()
	plt.plot(E_R * 1e6, F**2, color = color)
	plt.grid(linestyle='--', alpha=0.4)
	plt.xlabel(r'$E_R (keV)$', fontsize=14)# xunits=radians)
	plt.ylabel(r'$F^2$', fontsize=14)
	plt.title('Helm Form Factor')
	plt.text(20, 2e-9, r'$m_\chi = 1$ TeV',
         bbox=dict(facecolor='white', edgecolor='white', alpha=1), fontsize=10)
	plt.yscale('log')
	plt.xlim([0, 1e3])
	plt.ylim([1e-9, 2])
	plt.legend(['Xe', 'Ge', 'Ar', 'Na'])

#plt.show()
plt.savefig('/afs/ciemat.es/user/t/tomas/CIEMAT/Dark_Matter/programas/figures/FF_Helm(Ar).png')


