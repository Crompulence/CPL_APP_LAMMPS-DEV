import os
import sys
import errno
import pytest
import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt

# Import postproclib
sys.path.insert(0, "./pyDataView/")
try:
    import postproclib as ppl
except ImportError:
    cmd = "git clone https://github.com/edwardsmith999/pyDataView.git ./pyDataView"
    downloadout = sp.check_output(cmd, shell=True)
    sys.path.insert(0, "./pyDataView")
    import postproclib as ppl

# Add python scripts to path and import required classes
sys.path.append('../python_scripts/')
from LAMMPS_Input import LAMMPS_Input
from OpenFOAM_Input import OpenFOAM_Input
from DragForce import DragForce, Stokes, DiFelice, Ergun

# Run coupled simulation as subprocess
def run_coupled(run_bash_script='run.sh'):
    try:
        cmd = './' + run_bash_script
        p = sp.Popen(cmd, 
            stdout=sp.PIPE, 
            stderr=sp.STDOUT, 
            shell=True)
        while p.poll() is None:
            l = p.stdout.readline()
            print(l.rstrip())
        print(p.stdout.read())
    except:
        print('Error running bash run script' + run_bash_script + ' in base directory.')
        raise

# Read print data for the top particle on column
def read_print_data(print_file):
    data = np.loadtxt(print_file, skiprows=1)
    t = data[:,0]
    xy = data[:,1]

    return t, xy

# Extract the input parameters from DEM script for LAMMPS and OpenFOAM case
# directory. Check that they are consistent.
def get_input_parameters(md_input_file='./lammps/column.in', cfd_case_dir='./openfoam/'):
    # The porosity is hardwired to the value of a simple cubic array.
    epsf = 0.476401224402
    
    # LAMMPS script
    mObj = LAMMPS_Input(input_file=md_input_file)
    dragModel = mObj.read_variable('dragModel')
    muf = mObj.read_variable('dynamic_viscosity')
    rhof = mObj.read_variable('fluid_density')
    dp = mObj.read_variable('diameter')
    xlo = mObj.read_variable('xlo')
    xhi = mObj.read_variable('xhi')
    ylo = mObj.read_variable('ylo')
    yhi = mObj.read_variable('yhi')
    zlo = mObj.read_variable('zlo')
    zhi = mObj.read_variable('zhi')

    xyz_orig = np.array([xlo, ylo, zlo], order='F', dtype=np.float64)
    xyzL = np.array([xhi, yhi, zhi], order='F', dtype=np.float64)

    # OpenFOAM case 
    cObj = OpenFOAM_Input(case_dir=cfd_case_dir)
    _, Uf = cObj.read_0_file('./openfoam/0/Ub', 'inlet')
    Uf = Uf[1]
    rhof_ = cObj.read_constant_file('./openfoam/constant/transportProperties', 'rhob')
    nuf_ = cObj.read_constant_file('./openfoam/constant/transportProperties', 'nub')
    muf_ = nuf_*rhof_
    assert(rhof == rhof_)
    assert(muf == muf_)

    return dragModel, muf, rhof, Uf, dp, epsf, xyzL, xyz_orig

# Set the input parameters for the simulations. At present, only the inlet
# fluid velocity can be adjusted in the OpenFOAM input, while in the LAMMPS
# input, the drag model and movement type can be adjusted.
def set_input_parameters(Uf, dragModel, movementType, md_input_file='./lammps/column.in', cfd_case_dir='./openfoam/'):
    mObj = LAMMPS_Input(input_file=md_input_file)
    mObj.write_variable(variable_name='dragModel', variable_value=dragModel)
    mObj.write_variable(variable_name='movementType', variable_value=movementType)

    cObj = OpenFOAM_Input(case_dir=cfd_case_dir)
    cObj.write_0_file('./openfoam/0/Ub', 'inlet', Uf, isScalar=False, axis_val=1)

# Calculate the analytical solution for the pressure at the inlet, assuming
# that the pressure at the outlet is zero
def analytical_pressure(dragModel, muf, rhof, Uf, dp, epsf, xyzL, xyz_orig):
    # Obtain the pressure gradient across sample, and the inlet pressure
    # (assuming that the outlet pressure is zero)
    if dragModel == 'Drag' or dragModel == 'Stokes':
        fObj = Stokes(muf=muf, dp=dp)
    elif dragModel == 'DiFelice':
        fObj = DiFelice(muf=muf, rhof=rhof, Uf=Uf/epsf, dp=dp, Vp=0., epsf=epsf)
    elif dragModel == 'Ergun':
        fObj = Ergun(muf=muf, rhof=rhof, Uf=Uf/epsf, dp=dp, Vp=0., epsf=epsf)
    else:
        raise('Unknown drag force model specified')

    gradP = fObj.calculate_pressure_gradient(epsf=epsf, Uf=Uf/epsf, Vp=0., dp=dp)
    pSol = gradP*(xyzL[1] - xyz_orig[1])

    return pSol
    
# Obtain the pressure profile from the numerical simulation
def get_pressure_profile(CFD_folder='./openfoam/', pressure_var='p', axis_val=1):
    PPObj = ppl.OpenFOAM_PostProc(CFD_folder)
    pObj = PPObj.plotlist[pressure_var]
    h, p = pObj.profile(axis=axis_val, startrec=pObj.maxrec, endrec=pObj.maxrec)
    p = np.squeeze(p)

    return h, p

# Plot pressure profile obtained from numerical simulation and analytical
# solution. Save the file in the results directory (which is created if
# required) and also save the data used for plotting as .npz file.
def plot_pressure(h, p, pSol, xyz_orig, xyzL, file_name='fig_pressure'):
    plt.plot(h, p, 'r-')
    plt.plot([xyz_orig[1],xyzL[1]], [pSol,0.], 'k--')
    plt.xlabel('Height (cm)')
    plt.ylabel('Pressure (0.1Pa)')
    plt.legend(('Numerical', 'Analytical'))
    plt.tight_layout()
    if not os.path.exists(os.path.dirname(file_name)):
        try:
            os.makedirs(os.path.dirname(file_name))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
    plt.savefig(file_name + '.png')
    np.savez(file_name + '.npz', h, p, pSol, xyz_orig, xyzL)
    plt.close()

# Compare the pressure along the length of the sample at the end of the
# simulation for a specified relative error.
def compare_pressure(h, p, pSol, xyz_orig, xyzL, tol=0.01):
    pSolProfile = np.interp(h, [xyz_orig[1],xyzL[1]], [pSol,0.])
    for i in range(len(pSolProfile)):
        err = abs((pSolProfile[i] - p[i])/pSolProfile[i] <= tol)
        assert err, ('For h = {:.2f}, the obtained pressure of {:.6f} does not match analytical'.format(h[i], p[i])
                     +' solution {:.6f} within {:.2f}% relative error'.format(pSolProfile[i], tol*100))

# ----- Main ----- #
@pytest.fixture()
def setup():

    # Run coupled simulation
    run_coupled()

    # Extract input parameters
    dragModel, muf, rhof, Uf, dp, epsf, xyzL, xyz_orig = get_input_parameters()

    # Load print data
    # t, xy = read_print_data('lammps/print_column.txt')

    # Calculate analytical pressure profile
    pSol = analytical_pressure(dragModel, muf, rhof, Uf, dp, epsf, xyzL, xyz_orig)

    # Extract velocity, pressure and eps profile (at end of simulation only)
    h, p = get_pressure_profile()

    # Plot the results
    plot_pressure(h, p, pSol, xyz_orig, xyzL,
        file_name='./results/fig_pressure_Uf_{}_{}'.format(Uf, dragModel))

    return h, p, pSol, xyz_orig, xyzL

movementTypes = ['fixed', 'moving']
dragModels = ['DiFelice', 'Ergun']
Uf_values = [0.1, 0.2, 0.3, 0.4, 0.5]
@pytest.mark.parametrize('movementType', movementTypes)
@pytest.mark.parametrize('dragModel', dragModels)
@pytest.mark.parametrize('Uf', Uf_values)
def test_pressure(setup, Uf, dragModel, movementType):
    set_input_parameters(Uf, dragModel, movementType)
    h, p, pSol, xyz_orig, xyzL = setup
    compare_pressure(h, p, pSol, xyz_orig, xyzL, tol=0.02)

if __name__ == "__main__":
    h, p, pSol, ylo, yhi = setup()
