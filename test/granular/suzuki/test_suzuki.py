import os
import sys
import errno
import pytest
import subprocess as sp
import numpy as np
import time

# Add python scripts to path and import required classes
sys.path.append('../python_scripts/')
from LAMMPS_Input import LAMMPS_Input, LAMMPS_Writer
from DragForce import DiFelice, Ergun

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
def get_input_parameters(md_input_file='./lammps/suzuki.in'):    
    mObj = LAMMPS_Input(md_input_file)

    return mObj

# Set the input parameters for the simulations. At present, only the particle
# diameter and drag force model can be adjusted. Both these only apply to the
# LAMMPS input.
def set_input_parameters(Uf, dragModel, md_input_file='./lammps/suzuki.in'):
    LAMMPS_Writer(md_input_file, 'fluid_velocity', Uf)
    LAMMPS_Writer(md_input_file, 'dragModel', dragModel)

# Calculate the analytical displacement
def analytical_displacement(mObj):
    
    dragModel = mObj.dragModel
    muf = mObj.dynamic_viscosity
    rhof = mObj.fluid_density
    Uf = mObj.fluid_velocity
    epsf = mObj.porosity
    dp = mObj.diameter
    rhop = mObj.density
    kn = mObj.kn
    ylo = mObj.ylo
    yhi = mObj.yhi
    s = mObj.lattice_scale
    g = -mObj.gravity 

    # Obtain the pressure gradient across sample, and the inlet pressure
    # (assuming that the outlet pressure is zero)
    if dragModel == 'DiFelice':
        fObj = DiFelice(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=Uf/epsf, Vp=0.)
    elif dragModel == 'Ergun':
        fObj = Ergun(muf=muf, rhof=rhof, epsf=epsf, dp=dp, Uf=Uf/epsf, Vp=0.)
    else:
        raise ValueError('Unknown drag force model specified')

    fdi = fObj.calculate_drag_force(Uf=Uf/epsf, Vp=0.)
    volume = (np.pi/6)*dp**3
    fi = volume*rhop*g - volume*rhof*g + fdi/epsf
    H = yhi-ylo
    N = H/s
    disp = 0.5*N*(N+1)*fi/kn
    xySol = (yhi - 0.5*dp) + disp
    
    return xySol

def read_print_data(xy0, print_file='./lammps/print_suzuki.txt'):
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
    t = data[:,0]
    xy = data[:,1]

    # Append initial values
    t = np.insert(t, 0, 0)
    xy = np.insert(xy, 0, xy0)
    
    return t, xy

# Plot displacement and velocity profile obtained from numerical simulation
# and analytical solution. Save the file in the results directory (which is
# created if required) and also save the data used for plotting as .npz file.
def plot_displacement(t, xy, xySol, file_name='./fig'):
    # import matplotlib
    import matplotlib.pyplot as plt

    # Plot displacement
    plt.plot(t, xy, 'r-')
    plt.plot(t, np.ones_like(t)*xySol, 'k--')
    plt.xlabel('Time (s)')
    plt.ylabel('Position (cm)')
    plt.legend(('Numerical', 'Analytical'))
    plt.tight_layout()
    plt.ticklabel_format(useOffset=False)
    if not os.path.exists(os.path.dirname(file_name)):
        try:
            os.makedirs(os.path.dirname(file_name))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
    plt.savefig(file_name + '.png')
    np.savez(file_name + '.npz', t=t, xy=xy, xySol=xySol)
    plt.close()

# Test final displacement matches analytical solution.
def compare_displacement(xy, xySol, tol):
    err = (abs(xySol -xy[-1])/xySol) <= tol
    assert err, ('Final displacement of {:.6f} does not match analytical'.format(xy[-1])
                 +' solution of {:.6f} within {:.2f}% relative error.'.format(xySol_upward, tol*100))


# ----- Main ----- #
dragModels = ['DiFelice', 'Ergun']
Uf_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
@pytest.mark.parametrize('Uf', Uf_values)
@pytest.mark.parametrize('dragModel', dragModels)
def test_displacement(Uf, dragModel, plot_results=False):

    # Set input parameters
    set_input_parameters(Uf, dragModel)

    # Run coupled simulation
    run_coupled()
    if not os.path.exists('lammps/print_suzuki.txt'):
        print('Attempting to re-run coupled run.')
        run_coupled()

    # Extract input parameters
    mObj = get_input_parameters()

    # Load print data
    xy0 = mObj.yhi - 0.5*mObj.diameter
    t, xy = read_print_data(xy0)
    
    # Extract input parameters from lammps input script
    xySol = analytical_displacement(mObj)

    # Plot the results
    if plot_results:
        plot_displacement(t, xy, xySol,
            file_name='./results/fig_Uf_{}_{}'.format(Uf, dragModel))

    # Test analytical solution
    compare_displacement(xy, xySol, tol=0.02)
    
if __name__ == "__main__":
    test_displacement(Uf=0.2, dragModel='DiFelice', plot_results=True)
