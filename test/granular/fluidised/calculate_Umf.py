import sys
import numpy as np 
import matplotlib.pyplot as plt

# Add python scripts to path and import required classes
sys.path.append('../python_scripts/')
from LAMMPS_Input import LAMMPS_Input

# Extract relevant parameters from the LAMMPS input script.
mObj = LAMMPS_Input('./lammps/fluidised.in')
dp = mObj.diameter
epsf = mObj.porosity
rhop = mObj.density
rhof = mObj.fluid_density
muf = mObj.dynamic_viscosity
g = mObj.gravity

# Range of inlet (superficial) velocity
Umf = np.linspace(0.5, 2.0, 100);

# Effective weight term
Ar = np.ones_like(Umf)*(1 - epsf)*(rhop - rhof)*g

# Driving fluid term (Di Felice)
Re = rhof*dp*Umf/muf
chi = 3.7 - 0.65*np.exp(-0.5*((1.5 - np.log10(Re))**2))
Cd = (0.63 + 4.8/np.sqrt(Re))**2
gradP_DiFelice = (3*Cd*rhof/(4*dp))*(1 - epsf)*(epsf**(-1-chi))*(Umf**2)

# Driving fluid term (Ergun)
gradP_Ergun = (150*muf/(dp**2))*((1-epsf)**2/(epsf**3))*Umf + (1.75*rhof/dp)*((1-epsf)/(epsf**3))*(Umf**2)

# Create plot
plt.plot(Umf, gradP_DiFelice, 'b-')
plt.plot(Umf, gradP_Ergun, 'r-')
plt.plot(Umf, Ar, 'k--')
plt.xlabel('Superficial Inlet Velocity (cm/s)')
plt.ylabel('Pressure Gradient (Ba/cm)')
plt.legend(('Di Felice', 'Ergun'))
plt.tight_layout()
plt.ticklabel_format(useOffset=False)
plt.savefig('fluidisation_limit.png')
plt.close()