
import pytest
from cplpy import run_test, prepare_config
import subprocess as sp
import os
import glob
import numpy as np

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def get_subprocess_error(e):
    print("subprocess ERROR")
    import json
    error = json.loads(e[7:])
    print((error['code'], error['message']))


def runcmd(cmd):
    try:
        p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
        pcommout = p.communicate()
        output = pcommout[0].decode("utf-8")
        error = pcommout[1].decode("utf-8")
        if p.returncode != 0: 
            print("returncode", p.returncode)
            print("Stdout = ", output)
            print("Stderror = ", error)
        #run = sp.check_output(cmd, stderr=sp.STDOUT, shell=True).decode("utf-8")
    except sp.CalledProcessError as e:
        print("Stdout = ", e.stdout)
        print("Stderror = ", e.stderr)
        if e.output.startswith(b'error: {'):
            get_subprocess_error(e.output)
        raise

    return output+error


def execute(cmd, blocking=True):
    """
    Outputs results as they are generated, example usage:

    from __future__ import print_function

    for path in execute(["locate", "a"]):
        print(path, end="")

    """
    popen = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.STDOUT, shell=True, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line 
    popen.stdout.close()
    if blocking:
        return_code = popen.wait()

    if return_code:
        raise sp.CalledProcessError(return_code, cmd)


MD_EXEC = "lmp_cpl"
CFD_EXEC = "./CFD_single_ball.py"
TEST_DIR = os.path.dirname(os.path.realpath(__file__))

@pytest.fixture(scope="module")
def clean_dir():

    print("Cleaning directory")
    #Try to setup code
    with cd(TEST_DIR):
        try:
            clean = sp.check_output("rm -f " + "./thermo_output* " 
                                             + "./log.lammps* " 
                                             + "./debug.vels" 
                                             + " " + MD_EXEC.split("/")[-1], shell=True)
        except sp.CalledProcessError as e:
            if e.output.startswith(b'error: {'):
                get_subprocess_error(e.output)
            raise

    print(clean)

    return clean

@pytest.fixture(scope="module")
def build_case():

    print(("Building LAMMPS in ", TEST_DIR))
    #Try to setup code
    with cd(TEST_DIR):
        try:
            build = sp.check_output("./build.sh", shell=True)
        except sp.CalledProcessError as e:
            if e.output.startswith('error: {'):
                get_subprocess_error(e.output)
            raise

    return build


@pytest.fixture(scope="module")
def build_run(clean_dir, build_case):
    try:
        #clean = clean_dir()
        if os.path.isfile(MD_EXEC):
            print(("MD executable is ", MD_EXEC))
        #else:
        #    build = build_case()
    except sp.CalledProcessError:
        print("Build Failed")

def run_case(mdprocs):

    print(("Running case ", TEST_DIR))
    #Try to run code
    cmd = ('cplexec -p -Mv -m ' + str(mdprocs) + ' "' + MD_EXEC + ' -in single.in" ' + ' -c 1 ' +  CFD_EXEC)
    #cmd = ('mpiexec -n ' + str(mdprocs) + ' ' + MD_EXEC + ' -in ./single.in' + ' : -n 1 python ' +  CFD_EXEC)
    print(cmd)
    with cd(TEST_DIR):
        run = runcmd(cmd)

    return run

@pytest.mark.parametrize("mdprocs", [1, 2, 4])
def test_gravity(clean_dir, build_run, mdprocs):

    #Check vs analystical solution for gravity
    import bouncing
    run = run_case(mdprocs)

    with cd(TEST_DIR):
        error = bouncing.check_bouncing_error_vs_gravity()
        for e in error[0,1,:]:
            assert np.abs(e) < 1e-11

@pytest.mark.parametrize("mdprocs", [1, 2, 4])
def test_regression(clean_dir, build_run, mdprocs):

    #Check vs analystical solution for gravity
    import bouncing as b
    run = run_case(mdprocs)

    with cd(TEST_DIR):

        t, z, v, f = b.read_data(logfile='./log.lammps', 
                                 datafile='./thermo_output.txt')

        t_reg, z_reg, v_reg, f_reg = b.read_data(logfile='./regression_data/log.lammps', 
                                     datafile='./regression_data/thermo_output.txt')

        for i in range(z.shape[0]):
            assert np.abs((z[i]-z_reg[i])/z_reg[i]) < 1e-12
            assert np.abs((v[i]-v_reg[i])/z_reg[i]) < 1e-12
            assert np.abs((f[i]-f_reg[i])/z_reg[i]) < 1e-12

