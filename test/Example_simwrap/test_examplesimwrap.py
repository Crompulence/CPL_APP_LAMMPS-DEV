import pytest
import os
import sys
import numpy as np
import subprocess as sp

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)



sys.path.insert(0, "./SimWrapPy/")
try:
    import simwraplib as swl
except ImportError:
    cmd = "git clone https://github.com/edwardsmith999/SimWrapPy.git ./SimWrapPy"
    downloadout = sp.check_output(cmd, shell=True)
    print(downloadout)
    sys.path.insert(0, "./SimWrapPy")
    import simwraplib as swl

TEST_DIR = os.path.dirname(os.path.realpath(__file__))

gbase = 9.81
dtbase = 1e-5
params = []
for i in range(-10,5):
    ratio = 1. + 0.05*i
    params.append([dtbase/ratio, ratio*gbase])

@pytest.mark.parametrize("timestep, g", params)
def test_newtest(timestep, g):

    # Inputs that are the same for every thread
    basedir = TEST_DIR
    srcdir =  None
    executable = "/../../bin/lmp_cpl"
    inputfile = "/single.in"
    initstate = "/Bouncing_Ball.lj"
    rundir = TEST_DIR + "/run" # + "dt_" + str(np.round(timestep,2)) + "_g_" + str(np.round(g,2))

    #Clean previous result
    clean = sp.check_output("rm -f " + rundir + "./thermo_output* " 
                                     + rundir + "./log.lammps* " 
                                     + rundir + "./debug.vels" , shell=True)

    with cd(TEST_DIR):

        #Setup a LAMMPS run object
        lmps = swl.LammpsRun(None, basedir, rundir,
                             executable, inputfile, initstate=initstate,
                             inputchanges={"timestep":timestep})

        #Setup a mock script
        mockscript = "./CFD_single_ball.py"
        mock = swl.ScriptRun(rundir, mockscript, inputchanges={"g = ":g})

        #Setup a coupled run
        run = swl.CPLRun(None, basedir, rundir, [lmps, mock],
                         inputfile="/cpl/COUPLER.in")

        #Run the case
        run.setup()
        run.execute(blocking=True, print_output=False)

    #Check error
    import bouncing
    error = bouncing.check_bouncing_error_vs_gravity(g=g, logfile=rundir + '/log.lammps', 
                                                     datafile=rundir + '/thermo_output.txt')
    for e in error[0,1,:]:
        assert np.abs(e) < 1e-11
