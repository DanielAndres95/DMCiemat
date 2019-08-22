__author__ = ' AUTOR: Tomás Sánchez Sánchez-Pastor'

# Libraries:
import numpy as np
import math as math
import matplotlib.pyplot as plt
from scipy.special import erf

# Plot Settings:
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
colors = ['forestgreen', 'b', 'r', 'k']

# Constants:
mx_list = [10, 25, 50, 100] # GeV
Z = 18 
N = 22
A = Z + N
M = 931.5*A*1e-3 # GeV
rho = 0.3 # GeV cm⁻³
vo = 230 * 1e5 *3600 * 24 # cm day⁻¹
beta = vo / (1e2 * 3600 * 24 * 3e8)
No = 6.022e23*1e3 / A # kg⁻¹
vesc = 600 * 1e5 * 3600 *24 # cm day⁻¹
hc = 197.3e-3 # GeV fm
theta = np.linspace(0, math.pi, num=1e3, endpoint=False)
# Form Factor definitions:
c = 1.23*A**(1/3) - 0.60 # fm
s = 0.9 # fm
a = 0.52 # fm
sigma_n = 1e-45 # cm^-2

print("Nucleus with A = ", A)

# Load the data from DMTools:
data = np.loadtxt('Rec_Spec_100_5000.txt', delimiter=';').T

def kinematics(mx, M):
    mr = (mx*M/(mx + M)) # Masa reducida
    mrn = (mx*(M/A) / (mx + M/A)) # Masa reducida nucleon
    r = 4*mr**2/(mx*M)
    E_o = 1/2 * mx * beta**2
    E_R = 1/2 * E_o * r * (1 - np.cos(theta))
    q = np.sqrt(2 * M * E_R)
    output = {'mr': mr,
             'mrn': mrn,
             'Eo': E_o,
              'q': q,
              'r': r,
             'ER': E_R}
    return output

def form_factor(dic_kinematics):
    q = dic_kinematics['q']
    rn = np.sqrt(c**2 + 7/3*math.pi**2*a**2 - 5*s**2)
    j1 = (np.sin(q*rn/hc) - q*rn/hc*np.cos(q*rn/hc)) / (q*rn/hc)**2
    FF = 3 * j1 * np.exp(-(q*s/hc)**2/2) / (q*rn/hc)
    return FF

def compute_Ro(dic_kinematics):
    mr = dic_kinematics['mr']
    mrn = dic_kinematics['mrn']
    Ro = 2/math.sqrt(math.pi) * No*rho/mx * (mr/mrn)**2 * sigma_n * A * vo
    return Ro

def factors_of_dR(dic_kinematics):
    ER = dic_kinematics['ER']
    Eo = dic_kinematics['Eo']
    r = dic_kinematics['r']
    t = np.linspace(0, 366)
    
    k0_k1 = ( erf(vesc/vo) - 2/np.sqrt(math.pi) * vesc/vo * np.exp(- (vesc/vo)**2) )**(-1)
    v_E = ( 232 + 15*np.cos( 2*math.pi*(t - 152.5)/365.25 ) ) * 1e5 *3600 * 24 #cm day⁻¹
    v_min = np.sqrt( ER/(Eo*r) )*vo # cm day⁻¹
    
    output = {'k0_k1': k0_k1,
             'vE': np.max(v_E),
             'vmin': v_min}
    return output

def compute_dRdER(dic_kinematics, FF, Ro, factors_of_dR):
    r = dic_kinematics['r']
    Eo = dic_kinematics['Eo']
    k0_k1 = factors_of_dR['k0_k1']
    vE = factors_of_dR['vE']
    vmin = factors_of_dR['vmin']
    
    output = k0_k1 * Ro/( Eo*r ) * ( np.sqrt(math.pi)*vo/( 4*vE ) * ( erf( (vmin + vE)/vo ) - erf( (vmin - vE)/vo ) ) -
                                   np.exp( -( vesc/vo )**2 ) ) * FF**2 # events/kg/day/GeV
    return output

def Differential_Recoil_Spectrum(Wimp_mass, Nucleus_mass):
    kinematics_variables = kinematics(mx, M)
    FF = form_factor(kinematics_variables)
    Ro = compute_Ro(kinematics_variables)
    k0_k1 = ( erf(vesc/vo) - 2/np.sqrt(math.pi) * vesc/vo * np.exp(- (vesc/vo)**2) )**(-1)
    dR_factors = factors_of_dR(kinematics_variables)
    output = compute_dRdER(kinematics_variables, FF, Ro, dR_factors) * 1e-6
    print('dEdER units: events / kg / keV / day')
    return output, kinematics_variables['ER']

i = 0
fig, ax = plt.subplots(num=1, nrows=1, ncols=1, figsize=(8, 6))
for mx in mx_list:
    dRdER, ER = Differential_Recoil_Spectrum(mx, M)

    # Plots:
    cl = colors[i]
    i += 1
    plt.plot(ER*1e6, dRdER, color=cl, lw=2)
    plt.xlabel('$E_R\ (keV)$', fontsize=16)
    plt.ylabel('$Events / keV / kg / day$', fontsize=16)
    plt.xlim([0, 30])
    plt.yscale('log')
    plt.grid(alpha=.3)

plt.show()
