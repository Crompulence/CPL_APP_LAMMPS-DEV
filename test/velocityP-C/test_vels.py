#!/usr/bin/env python2
import pytest
import numpy as np
from cplpy import run_test, prepare_config, get_test_dir
import os
import re
from testutils import check_vels

# -----Velocities TESTS-----

# EXPLANATION:

MD_FNAME = "lammps_vels.in"
MD_ARGS = "-in " + MD_FNAME
MD_EXEC = "lmp_cpl"
CFD_FNAME = "dummyCFD_vels.py"
CFD_ARGS = CFD_FNAME
CFD_EXEC = "python2"
TEST_TEMPLATE_DIR = os.path.join(get_test_dir(), "templates")
TEST_DIR = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture()
def prepare_config_fix(tmpdir):
    prepare_config(tmpdir, TEST_DIR, MD_FNAME, CFD_FNAME)

# -----VELOCITY TESTS-----

# EXPLANATION: See README-test located in this folder.


@pytest.mark.parametrize("cfdprocs, mdprocs, err_msg", [
                         ((3, 3, 3), (3, 3, 3),  ""),
                         ((1, 1, 1), (3, 3, 3),  "")])
def test_velocitiesP2C(prepare_config_fix, cfdprocs, mdprocs, err_msg):
    MD_PARAMS = {"lx": 300.0, "ly": 300.0, "lz": 300.0}
    MD_PARAMS["npx"], MD_PARAMS["npy"], MD_PARAMS["npz"] = mdprocs

    CFD_PARAMS = {"lx": 300.0, "ly": 300.0, "lz": 300.0,
                  "ncx": 15, "ncy": 15, "ncz": 15, }
    CFD_PARAMS["npx"], CFD_PARAMS["npy"], CFD_PARAMS["npz"] = cfdprocs

    CONFIG_PARAMS = {"cfd_bcx": 1, "cfd_bcy": 1, "cfd_bcz": 1,
                     "olap_xlo": 1, "olap_xhi": 15,
                     "olap_ylo": 1, "olap_yhi": 5,
                     "olap_zlo": 1, "olap_zhi": 15,
                     "cnst_xlo": 1, "cnst_xhi": 15,
                     "cnst_ylo": 5, "cnst_yhi": 5,
                     "cnst_zlo": 1, "cnst_zhi": 15,
                     "bndry_xlo": 1, "bndry_xhi": 15,
                     "bndry_ylo": 1, "bndry_yhi": 1,
                     "bndry_zlo": 1, "bndry_zhi": 15,
                     "tstep_ratio": 5, }

    correct = run_test(TEST_TEMPLATE_DIR, CONFIG_PARAMS, MD_EXEC, MD_FNAME, MD_ARGS,
                       CFD_EXEC, CFD_FNAME, CFD_ARGS, MD_PARAMS, CFD_PARAMS, err_msg, True)
    if correct:
        check_vels(1e-6, steps=2)
