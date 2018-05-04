import pytest
import os
import sys
import numpy as np
import subprocess as sp

sys.path.append("/home/es205/codes/python/SimWrapPy/")
try:
    import simwraplib as swl
except ImportError:
    print("Downloading simwraplib")
    import subprocess as sp
    downloadout = sp.check_output("git clone https://github.com/edwardsmith999/SimWrapPy.git ./SimWrapPy", shell=True)
    print(downloadout)
    sys.path.append("./SimWrapPy/")
    import simwraplib as swl

TEST_DIR = os.path.dirname(os.path.realpath(__file__))

params = list(np.arange(1, 3)*1e-5)

@pytest.mark.parametrize("timestep", params)
def test_newtest(timestep):

    # Inputs that are the same for every thread
    basedir = TEST_DIR
    srcdir =  None
    executable = "/../../bin/lmp_cpl"
    inputfile = "/single.in"
    initstate = "/Bouncing_Ball.lj"
    rundir = TEST_DIR + "/run"

    #Clean previous result
    clean = sp.check_output("rm -f " + rundir + "./thermo_output* " 
                                     + rundir + "./log.lammps* " 
                                     + rundir + "./debug.vels" , shell=True)

    #Setup a LAMMPS run object
    lmps = swl.LammpsRun(None, basedir, rundir,
                         executable, inputfile, initstate=initstate,
                         inputchanges={"timestep":timestep})

    #Setup a mock script
    mockscript = "./CFD_single_ball.py"
    mock = swl.ScriptRun(rundir, mockscript)

    #Setup a coupled run
    run = swl.CPLRun(None, basedir, rundir, [lmps, mock],
                     inputfile="/cpl/COUPLER.in")

    #Run the case
    run.setup()
    run.execute(blocking=True)

    #Check error
    import bouncing
    error = bouncing.check_bouncing_error_vs_gravity(logfile=rundir + '/log.lammps', 
                                                     datafile=rundir + '/thermo_output.txt')
    for e in error[0,1,:]:
        assert np.abs(e) < 1e-11
