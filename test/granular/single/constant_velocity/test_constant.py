import os
import sys
import errno
import pytest
import subprocess as sp
import numpy as np
import time

# Add python scripts to path and import required classes
sys.path.append('../../python_scripts/')
from LAMMPS_Input import LAMMPS_Input, LAMMPS_Writer
from MOCK_Input import MOCK_Input
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
        raise RuntimeError('Error running bash run script' + run_bash_script + ' in base directory.')
    p.wait()

# Extract the input parameters from DEM script for LAMMPS and OpenFOAM case
# directory. Check that they are consistent.
def get_input_parameters(md_input_file='./lammps/constant.in', cfd_input_file='./CFD_dummy_constant.py'):    
    mObj = LAMMPS_Input(md_input_file)
    cObj = MOCK_Input(cfd_input_file)

    # Calculate porosity and add to mObj
    Vc = ((cObj.xyzL[0]-cObj.xyz_orig[0])/cObj.ncxyz[0])*((cObj.xyzL[1]-cObj.xyz_orig[1])/cObj.ncxyz[1])*((cObj.xyzL[2]-cObj.xyz_orig[2])/cObj.ncxyz[2])
    mObj.epsf = (Vc - (np.pi/6)*(mObj.diameter**3))/Vc

    return mObj

# Set the input parameters for the simulations. At present, only the particle
# diameter and drag force model can be adjusted. Both these only apply to the
# LAMMPS input.
def set_input_parameters(dp, dragModel, vy0, md_input_file='./lammps/constant.in'):
    LAMMPS_Writer(md_input_file, 'diameter', dp)
    LAMMPS_Writer(md_input_file, 'dragModel', dragModel)
    LAMMPS_Writer(md_input_file, 'vy0', vy0)

# Calculate the analytical displacement and velocity profile of the particle,
# along with the terminal velocity. This is only applicable for the Stokes
# drag law, at present.
def analytical_force(mObj):
    
    dragModel = mObj.dragModel
    muf = mObj.dynamic_viscosity
    rhof = mObj.fluid_density
    dp = mObj.diameter
    vy0 = mObj.vy0
    epsf = mObj.epsf
    
    if dragModel == 'Drag' or dragModel == 'Stokes':
        fObj = Stokes(muf=muf, epsf=epsf, dp=dp)
    elif dragModel == 'DiFelice':
        fObj = DiFelice(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=0., Vp=vy0)
    elif dragModel == 'Ergun':
        fObj = Ergun(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=0., Vp=vy0)
    else:
        raise ValueError('Unknown drag force model specified')

    fySol = fObj.calculate_drag_force(Uf=0., Vp=vy0)

    return fySol

# Read print data for the top particle on column
def read_print_data(print_file='./lammps/print_constant.txt'):
    # Try reading the print file. StopIteration error occurs with 'CPL_init
    # has been called more than once. Returning same COMM' error during
    # coupled run causing the print_file to exist but be empty (and hence the
    # skip_header error). Temporary work around is to re-run the coupled
    # simulation after waiting for 3 seconds. Only attempt this re-run three
    # times.
    for i in range(3):
        try:
            data = np.genfromtxt(print_file, skip_header=1)
            break
        except StopIteration:
            print("genfromtxt read error, waiting 3 secs and will try again")
            time.sleep(3)
            run_coupled()
    
    # Extract data
    fy = data[:,3]

    return fy

# Compare the displacement profile with time for a specified relative error.
def compare_force(fy, fySol, tol=0.01):
    for i in range(len(fy)-1):
        err = abs((fySol - fy[i])/fySol) <= tol
        assert err, ('Imposed drag force of {:.6f} does not match analytical'.format(fy[i])
                + ' solution of {:.6f} within {:.2f}% relative error.'.format(fySol, tol*100))

# ----- Main ----- #
# dragModels = ['Drag', 'Stokes', 'DiFelice', 'Ergun']
dragModels = ['DiFelice']
dp_values = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10]
vy0_values = [0.5, 1.0, 1.5, 2.0, 2.5]
@pytest.mark.parametrize('vy0', vy0_values)
@pytest.mark.parametrize('dragModel', dragModels)
@pytest.mark.parametrize('dp', dp_values)
def test_force(dp, dragModel, vy0, plot_results=False):
    
    # Set input parameters
    set_input_parameters(dp, dragModel, vy0)
    
    # Run coupled simulation
    run_coupled()

    # Extract input parameters
    mObj = get_input_parameters()

    # Extract input parameters from lammps input script
    fySol = analytical_force(mObj)
    
    # Load print data
    fy = read_print_data()
    
    # Test analytical solution
    compare_force(fy, fySol, tol=0.02)

    # Save results
    if plot_results:
        file_name = './results/data_{}_dp_{}_vy0_{}.npz'.format(dragModel, dp, vy0)
        if not os.path.exists(os.path.dirname(file_name)):
            try:
                os.makedirs(os.path.dirname(file_name))
            except OSError as exc:
                if exc.errno != errno.EEXIST:
                    raise 
        np.savez(file_name, fy=fy[-2], fySol=fySol)

    # time.sleep(5)

if __name__ == "__main__":
    test_force(dp=0.10, dragModel='Ergun', vy0=0.5)
