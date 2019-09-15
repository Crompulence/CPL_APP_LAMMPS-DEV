import sys
import numpy as np
from mpi4py import MPI
from cplpy import CPL

# Add python scripts to path and import required classes
sys.path.append('../../python_scripts/')
from LAMMPS_Input import LAMMPS_Input
from DragForce import DragForce, Stokes, DiFelice, Ergun

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
mObj = LAMMPS_Input('./lammps/constant.in')
dragModel = mObj.dragModel
muf = mObj.dynamic_viscosity
rhof = mObj.fluid_density
dp = mObj.diameter
vy0 = mObj.vy0

CELL_VOLUME = ((xyzL[0]-xyz_orig[0])/ncxyz[0])*((xyzL[1]-xyz_orig[1])/ncxyz[1])*((xyzL[2]-xyz_orig[2])/ncxyz[2])
epsf = (CELL_VOLUME - (np.pi/6)*(dp**3))/CELL_VOLUME

if dragModel == 'Drag' or dragModel == 'Stokes':
    fObj = Stokes(muf=muf, epsf=epsf, dp=dp)
elif dragModel == 'DiFelice':
    fObj = DiFelice(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=0., Vp=vy0)
elif dragModel == 'Ergun':
    fObj = Ergun(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=0., Vp=vy0)
else:
    raise ValueError('Unknown drag force model specified')

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
for time in range(20):

    # print(time)

    # Send data: Zero send buffer. Set only the fluid velocity and gradP in the
    # y-direction for the entire column.
    # [Ux, Uy, Uz, gradPx, gradPy, gradPz, divTaux, divTauy, divTauz]
    sendbuf[:,:,:,:] = 0
    sendbuf[4,:,:,:] = -fObj.calculate_pressure_gradient(Uf=0., Vp=vy0)

    CPL.send(sendbuf)

    # Recieve data:
    # [UxSum, UySum, UzSum, FxSum, FySum, Fz6Sum, CdSum, VolSum] 
    recvbuf, ierr = CPL.recv(recvbuf)

CPL.finalize()
MPI.Finalize()
