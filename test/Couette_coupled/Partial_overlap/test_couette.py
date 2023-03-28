import pytest
import os
import sys
import numpy as np
import subprocess as sp
from check_LAMMPS_vs_analytical import check_LAMMPS_vs_Analytical


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

# Import symwraplib
sys.path.insert(0, "./SimWrapPy/")
try:
    import simwraplib as swl
except ImportError:
    cmd = "git clone https://github.com/edwardsmith999/SimWrapPy.git ./SimWrapPy"
    downloadout = sp.check_output(cmd, shell=True)
    print(downloadout)
    sys.path.insert(0, "./SimWrapPy")
    import simwraplib as swl

#Define test directory based on script file
TEST_DIR = os.path.dirname(os.path.realpath(__file__))

#Parameterise range of cases
params = []
#Uwall = [0.6, 0.7, 0.8, 0.9, 1.0]
Uwall = [0.6,0.7,0.8,0.9,1.0]
for u in Uwall:
    params.append(u)
@pytest.mark.parametrize("wallvel", params)
def test_newtest(wallvel):

    # Inputs that are the same for every thread
    basedir = TEST_DIR
    srcdir = None
    inputfile = "/one_wall_imaginary_sliding_coupled.in"
    rundir = TEST_DIR + "/run" + str(wallvel)
    executable = "lmp_cpl"

    with cd(TEST_DIR):

        #Setup a LAMMPS run object
        lmps = swl.LammpsRun(None, basedir, rundir,
                             executable, inputfile)

        #Setup a mock script
        mockscript = "./python_dummy/CFD_test_vs_Couette.py"
        mock = swl.ScriptRun(rundir, mockscript, inputchanges={"Uwall = ": wallvel})

        #Setup a coupled run
        run = swl.CPLRun(None, basedir, rundir, [lmps, mock],
                         inputfile="/cpl/COUPLER.in")

        #Run the case
        run.setup()
        run.execute(blocking=True, print_output=True, out_to_file=False, extra_cmds="-M -p")

        #Check results are correct
        check_LAMMPS_vs_Analytical(rundir, uwall=wallvel, tol=0.12)

            
