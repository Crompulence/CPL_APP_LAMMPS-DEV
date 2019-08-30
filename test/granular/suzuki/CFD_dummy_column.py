import sys
import numpy as np
from mpi4py import MPI
from cplpy import CPL

# Add python scripts to path and import required classes
sys.path.append('../python_scripts/')
from DragForce import DragForce, Stokes, DiFelice

# --------------------------------------------------------
# USER INPUT 
# --------------------------------------------------------
# Parameters of the cpu topology (cartesian grid)
npxyz = np.array([1, 1, 1], order='F', dtype=np.int32)
xyzL = np.array([0.1, 10.0, 0.1], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)
ncxyz = np.array([1, 100, 1], order='F', dtype=np.int32)

# Parameters for model (ensure that they are the same as lammps/column.in)
dragModel = 'DiFelice'
muf = 1.e-2
rhof = 1.0
Uf = 0.1
epsf = 0.4764
dp = 0.1
g = 981.
# --------------------------------------------------------
# --------------------------------------------------------

# initialise MPI
comm = MPI.COMM_WORLD

# Initialise CPL
CPL = CPL()
CFD_COMM = CPL.init(CPL.CFD_REALM)
cart_comm = CFD_COMM.Create_cart([npxyz[0], npxyz[1], npxyz[2]])
CPL.setup_cfd(cart_comm, xyzL, xyz_orig, ncxyz)

# Setup send and recv buffers
recvbuf, sendbuf = CPL.get_arrays(recv_size=8, send_size=9)

# Main time loop
for time in range(40000):

	# print(time)
	
	if time == 0:
		Vp = np.array([0.0])

	# Setup drag force object
	if dragModel == 'Drag' or dragModel == 'Stokes':
	    fObj = Stokes(muf=muf, dp=dp)
	elif dragModel == 'DiFelice':
	    if len(Vp) == 1:
		    fObj = DiFelice(muf=muf, rhof=rhof, Uf=Uf/epsf, dp=dp, Vp=Vp, epsf=epsf)
		    gradPy = fObj.calculate_pressure_gradient(epsf=epsf, Uf=Uf, Vp=Vp, dp=dp) + rhof*g
	    elif len(Vp) > 1:
	    	gradPy = np.zeros_like(Vp)
	    	for i in range(len(Vp)):
	    		fObj = DiFelice(muf=muf, rhof=rhof, Uf=Uf/epsf, dp=dp, Vp=Vp[i], epsf=epsf)
	    		gradPy[i] = fObj.calculate_pressure_gradient(epsf=epsf, Uf=Uf, Vp=Vp[i], dp=dp) + rhof*g
    		gradPy = np.reshape(gradPy, (1,100,1))
	    else:
	    	raise('Error with calculation of particle velocity')
	else:
	    raise('Unknown drag force model specified')

	# Send data: Zero send buffer. Set only the fluid velocity and gradP in the
	# y-direction for the entire column.
	# [Ux, Uy, Uz, gradPx, gradPy, gradPz, divTaux, divTauy, divTauz]
	sendbuf[:,:,:,:] = 0
	sendbuf[1,:,:,:] = Uf/epsf
	sendbuf[4,:,:,:] = -gradPy

	CPL.send(sendbuf)

	# Recieve data:
	# [UxSum, UySum, UzSum, FxSum, FySum, FzSum, CdSum, VolSum] 
	recvbuf, ierr = CPL.recv(recvbuf)

	# Calculate particle velocity field which will affect the obtained pressure gradient
	Vp = np.squeeze(recvbuf[1,:,:,:])/((np.pi/6)*(dp**3))
	# if time > 0 and (time % 100) == 0:
	# 	print('Output particle velocity for top and bottom particle')
	# 	print(Vp[1])
	# 	print(Vp[-1])

CPL.finalize()
MPI.Finalize()