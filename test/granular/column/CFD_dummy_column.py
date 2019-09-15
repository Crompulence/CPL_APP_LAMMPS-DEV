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
# Parameters of the cpu topology (cartesian grid)
npxyz = np.array([1, 1, 1], order='F', dtype=np.int32)
xyzL = np.array([0.1, 10.0, 0.1], order='F', dtype=np.float64)
xyz_orig = np.array([0.0, 0.0, 0.0], order='F', dtype=np.float64)
ncxyz = np.array([1, 100, 1], order='F', dtype=np.int32)
# --------------------------------------------------------
# USER INPUT - END
# --------------------------------------------------------

# Extract relevant parameters from the LAMMPS input script.
mObj = LAMMPS_Input(input_file='./lammps/resting.in')
dragModel = mObj.read_variable('dragModel')
muf = mObj.read_variable('dynamic_viscosity')
rhof = mObj.read_variable('fluid_density')
Uf = mObj.read_variable('fluid_velocity')
dp = mObj.read_variable('diameter')
rhop = mObj.read_variable('density')
g = mObj.read_variable('gravity')
g = -g
epsf = (0.1**3 - (np.pi/6)*dp**3)/(0.1**3)

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
    if time == 0:
        Vp = 0.


    print(Vp)

    if dragModel == 'Drag' or dragModel == 'Stokes':
        fObj = Stokes(muf=muf, dp=dp, epsf=epsf)
    elif dragModel == 'DiFelice':
        fObj = DiFelice(muf=muf, rhof=rhof, Uf=Uf/epsf, dp=dp, Vp=Vp, epsf=epsf)
    elif dragModel == 'Ergun':
        fObj = Ergun(muf=muf, rhof=rhof, Uf=Uf/epsf, dp=dp, Vp=Vp, epsf=epsf)
    else:
        raise('Unknown drag force model specified')

    # Send data: Zero send buffer. Set only the fluid velocity and gradP in the
    # y-direction for the entire column.
    # [Ux, Uy, Uz, gradPx, gradPy, gradPz, divTaux, divTauy, divTauz]
    sendbuf[:,:,:,:] = 0
    sendbuf[1,:,:,:] = Uf/epsf
    sendbuf[4,:,:,:] = rhof*g - fObj.calculate_pressure_gradient(epsf=epsf, Uf=Uf/epsf, Vp=Vp, dp=dp)

    CPL.send(sendbuf)

    # Recieve data:
    # [UxSum, UySum, UzSum, FxSum, FySum, FzSum, CdSum, VolSum] 
    recvbuf, ierr = CPL.recv(recvbuf)
    Vp = np.squeeze(recvbuf[1,:,:,:])/((np.pi/6)*(dp**3))
    Vp = Vp[np.nonzero(Vp)]
    if len(Vp) == 1:
        Vp = Vp[0]
    else:
        print('Multiple cells detected with particle velocity, even though it is single particle.')
        print(Vp)
        Vp = 0.

CPL.finalize()
MPI.Finalize()
