import sys
import numpy as np
from mpi4py import MPI
from cplpy import CPL

# Add python scripts to path and import required classes
sys.path.append('../python_scripts/')
from LAMMPS_Input import LAMMPS_Input
from DragForce import DiFelice, Ergun

# --------------------------------------------------------
# USER INPUT 
# --------------------------------------------------------
npxyz = [1, 1, 1]
xyzL = [0.1, 10.0, 0.1]
xyz_orig = [0.0, 0.0, 0.0]
ncxyz = [1, 100, 1]
# --------------------------------------------------------
# --------------------------------------------------------

# Parameters of the cpu topology (cartesian grid)
npxyz = np.array(npxyz, order='F', dtype=np.int32)
xyzL = np.array(xyzL, order='F', dtype=np.float64)
xyz_orig = np.array(xyz_orig, order='F', dtype=np.float64)
ncxyz = np.array(ncxyz, order='F', dtype=np.int32)

# Extract relevant parameters from the LAMMPS input script.
mObj = LAMMPS_Input('./lammps/suzuki.in')
dragModel = mObj.dragModel
muf = mObj.dynamic_viscosity
rhof = mObj.fluid_density
Uf = mObj.fluid_velocity
epsf = mObj.porosity
dp = mObj.diameter
g = -mObj.gravity

# initialise MPI
comm = MPI.COMM_WORLD

# Initialise CPL
CPL = CPL()
CFD_COMM = CPL.init(CPL.CFD_REALM)
cart_comm = CFD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]])
CPL.setup_cfd(cart_comm, xyzL, xyz_orig, ncxyz)

# Setup send and recv buffers
recvbuf, sendbuf = CPL.get_arrays(recv_size=8, send_size=9)

# Setup drag force object
if dragModel == 'DiFelice':
    fObj = DiFelice(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=Uf/epsf, Vp=0.)
elif dragModel == 'Ergun':
    fObj = Ergun(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=Uf/epsf, Vp=0.)
else:
    raise ValueError('Unknown drag force model specified')

# Main time loop
for time in range(400):

	# print(time)
	
	if time == 0:
		Vp = np.array([0.0])

	# Update drag coefficients for new particle velocities and calculate gradient
	# due to drag forces.
	if len(Vp) == 1:
		fObj.update_drag_coefficient(Uf=Uf/epsf, Vp=Vp)
		gradP = fObj.calculate_pressure_gradient(Uf=Uf/epsf, Vp=Vp)
	elif len(Vp) > 1:
		gradP = np.zeros_like(Vp)
		for i in range(len(Vp)):
			fObj.update_drag_coefficient(Uf=Uf/epsf, Vp=Vp[i])
			gradP[i] = fObj.calculate_pressure_gradient(Uf=Uf/epsf, Vp=Vp[i])
		gradP = np.reshape(gradP, (1,len(Vp),1))
	else:
		raise ValueError('Error with calculation of particle velocity')

	# Send data: Zero send buffer. Set only the fluid velocity and gradP in the
	# y-direction for the entire column.
	# [Ux, Uy, Uz, gradPx, gradPy, gradPz, divTaux, divTauy, divTauz]
	sendbuf[:,:,:,:] = 0
	sendbuf[1,:,:,:] = Uf/epsf
	sendbuf[4,:,:,:] = rhof*g - gradP

	CPL.send(sendbuf)

	# Recieve data:
	# [UxSum, UySum, UzSum, FxSum, FySum, FzSum, CdSum, VolSum] 
	recvbuf, ierr = CPL.recv(recvbuf)

	# Calculate particle velocity field which will affect the obtained pressure gradient
	Vp = np.squeeze(recvbuf[1,:,:,:])/((np.pi/6)*(dp**3))

CPL.finalize()
MPI.Finalize()
