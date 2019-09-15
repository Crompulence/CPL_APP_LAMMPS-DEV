import sys
import numpy as np
from mpi4py import MPI
from cplpy import CPL

# Add python scripts to path and import required classes
sys.path.append('../../python_scripts/')
from LAMMPS_Input import LAMMPS_Input
from DragForce import DragForce, Stokes

# --------------------------------------------------------
# USER INPUT - START 
# --------------------------------------------------------
npxyz = [1, 1, 1]
xyzL = [0.1, 10.0, 0.1]
xyz_orig = [0.0, 0.0, 0.0]
ncxyz = [1, 100, 1]
# --------------------------------------------------------
# USER INPUT - END
# --------------------------------------------------------

# Parameters of the cpu topology (cartesian grid)
npxyz = np.array(npxyz, order='F', dtype=np.int32)
xyzL = np.array(xyzL, order='F', dtype=np.float64)
xyz_orig = np.array(xyz_orig, order='F', dtype=np.float64)
ncxyz = np.array(ncxyz, order='F', dtype=np.int32)

# Extract relevant parameters from the LAMMPS input script.
mObj = LAMMPS_Input(input_file='./lammps/terminal.in')
muf = mObj.dynamic_viscosity
rhof = mObj.fluid_density
dp = mObj.diameter
rhop = mObj.density
g = -mObj.gravity

Vc = ((xyzL[0]-xyz_orig[0])/ncxyz[0])*((xyzL[1]-xyz_orig[1])/ncxyz[1])*((xyzL[2]-xyz_orig[2])/ncxyz[2])
epsf = (Vc - (np.pi/6)*(dp**3))/Vc

# Initialise drag force object
fObj = Stokes(muf=muf, dp=dp, epsf=epsf)

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
for time in range(400):

    # print(time)
    if time == 0:
        Vp = 0.

    # Send data: Zero send buffer. Set only the fluid velocity and gradP in the
    # y-direction for the entire column.
    # [Ux, Uy, Uz, gradPx, gradPy, gradPz, divTaux, divTauy, divTauz]
    sendbuf[:,:,:,:] = 0
    sendbuf[4,:,:,:] = rhof*g - fObj.calculate_pressure_gradient(Uf=0., Vp=Vp)

    CPL.send(sendbuf)

    # Recieve data:
    # [UxSum, UySum, UzSum, FxSum, FySum, FzSum, CdSum, VolSum] 
    recvbuf, ierr = CPL.recv(recvbuf)

    # Extract particle velocity for UySum = volume*Vp. Note that when
    # particles are detected in multiple cells (unsure why this occurs when no
    # overlap is considered), sum the value to obtain the single particle
    # velocity.
    Vp = np.squeeze(recvbuf[1,:,:,:])/((np.pi/6)*(dp**3))
    Vp = Vp[np.nonzero(Vp)]
    if len(Vp) == 1:
        Vp = Vp[0]
    elif len(Vp) == 0:
        print('No cells detected with particle velocity. Setting to zero.')
        Vp = 0.
    else:
        print('Multiple cells detected with particle velocity. Setting to sum of value.')
        Vp = np.sum(Vp)

CPL.finalize()
MPI.Finalize()
