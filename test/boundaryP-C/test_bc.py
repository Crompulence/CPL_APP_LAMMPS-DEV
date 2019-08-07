#!/usr/bin/env python2
import pytest
import numpy as np
from cplpy import run_test, prepare_config, get_test_dir, parametrize_file
import os, sys
import re

# -----Velocities TESTS-----

# EXPLANATION:

MD_FNAME = "lammps_bc.in"
MD_ARGS = "-in " + MD_FNAME
MD_EXEC = "lmp_cpl"
CFD_FNAME = "dummyCFD_bc.py"
CFD_ARGS = CFD_FNAME
CFD_EXEC = "python2"
TEST_TEMPLATE_DIR = os.path.join(get_test_dir(), "templates")
TEST_DIR = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture()
def prepare_config_fix(tmpdir):
    prepare_config(tmpdir, TEST_DIR, [MD_FNAME, CFD_FNAME, "config.cpl"])

# -----VELOCITY TESTS-----

# EXPLANATION: See README-test located in this folder.


@pytest.mark.parametrize("cfdprocs, mdprocs, ncells, ncy_olap, bin_dimension, err_msg", [
                         ((2, 1, 2), (2, 1, 2), (10, 15, 10),  5, 3, ""),
                         ((2, 1, 1), (2, 1, 2), (20, 15, 10),  5, 3, ""),
                         ((1, 1, 2), (1, 1, 2), (5, 15, 10), 3, 3, ""),
                         ((1, 3, 1), (1, 3, 1), (15, 15, 15), 5, 3, ""),
                         ((2, 1, 2), (2, 1, 2), (10, 15, 10),  2, 3, ""),
                         ((2, 1, 1), (2, 1, 2), (20, 15, 10),  2, 3, ""),
                         ((1, 1, 2), (1, 1, 2), (5, 15, 10), 3, 3, ""),
                         ((1, 3, 1), (1, 3, 1), (15, 15, 15), 2, 3, ""),
                         ((2, 1, 2), (2, 1, 2), (10, 15, 10),  5, 1, ""),
                         ((2, 1, 1), (2, 1, 2), (20, 15, 10),  5, 1, ""),
                         ((1, 1, 2), (1, 1, 2), (5, 15, 10), 3, 1, ""),
                         ((1, 3, 1), (1, 3, 1), (15, 15, 15), 5, 1, ""),
                         ((2, 1, 2), (2, 1, 2), (10, 15, 10),  2, 1, ""),
                         ((2, 1, 1), (2, 1, 2), (20, 15, 10),  2, 1, ""),
                         ((1, 1, 2), (1, 1, 2), (5, 15, 10), 3, 1, ""),
                         ((1, 3, 1), (1, 3, 1), (15, 15, 15), 2, 1, ""),
                         ((1, 1, 1), (1, 1, 1), (15, 15, 15), 5, 1, "")])
def test_velocitiesP2C(prepare_config_fix, cfdprocs, mdprocs, ncells, ncy_olap, bin_dimension, err_msg):
    MD_PARAMS = {"lx": 300.0, "ly": 300.0, "lz": 300.0, "bin_dimension": bin_dimension,
                  "ncx": ncells[0], "ncy": ncells[1], "ncz": ncells[2], }
    MD_PARAMS["npx"], MD_PARAMS["npy"], MD_PARAMS["npz"] = mdprocs

    CFD_PARAMS = {"lx": 300.0, "ly": 300.0, "lz": 300.0,
                  "ncx": ncells[0], "ncy": ncells[1], "ncz": ncells[2], }
    CFD_PARAMS["npx"], CFD_PARAMS["npy"], CFD_PARAMS["npz"] = cfdprocs

    CONFIG_PARAMS = {"cfd_bcx": 1, "cfd_bcy": 1, "cfd_bcz": 1,
                     "olap_xlo": 1, "olap_xhi": ncells[0],
                     "olap_ylo": 1, "olap_yhi": ncy_olap,
                     "olap_zlo": 1, "olap_zhi": ncells[2],
                     "cnst_xlo": 1, "cnst_xhi": ncells[0],
                     "cnst_ylo": ncy_olap, "cnst_yhi": ncy_olap,
                     "cnst_zlo": 1, "cnst_zhi": ncells[2],
                     "bndry_xlo": 1, "bndry_xhi": ncells[0],
                     "bndry_ylo": 1, "bndry_yhi": 1,
                     "bndry_zlo": 1, "bndry_zhi": ncells[2],
                     "tstep_ratio": 5, }
    MD_PARAMS.update(CONFIG_PARAMS)
    parametrize_file("lammps_bc.in", "lammps_bc.in", MD_PARAMS)
    parametrize_file("config.cpl", "config.cpl", MD_PARAMS)
    correct = run_test(TEST_TEMPLATE_DIR, CONFIG_PARAMS, MD_EXEC, MD_FNAME, MD_ARGS,
                       CFD_EXEC, CFD_FNAME, CFD_ARGS, MD_PARAMS, CFD_PARAMS, err_msg, True)
    if correct:
        compare_bc(1e-6, steps=3, bin_dimension=bin_dimension)


def name_with_step(fname, step):
    split_fname = re.split('(\W)', fname)
    dot_pos = split_fname.index(".")
    split_fname.insert(dot_pos , str(step))
    return "".join(split_fname)

def compare_bc(tol, lammps_fname="lammps_bc.dat", cfd_fname="cfd_bc.dat", steps=1, bin_dimension=3, mode="test"):

    assert bin_dimension == 1 or bin_dimension == 3
    # Line format of CFD script file -- > x y z vx vy vz
    for s in xrange(steps):
        with open(name_with_step(cfd_fname, s), "r") as cfd_file:
            # print cfd_fname
            cfd_lines = cfd_file.readlines()
        cfd_lines = [l[:-1].split(" ") for l in cfd_lines]
        cfd_cells = {}
        for l in cfd_lines:
            cfd_cells[(float(l[0]), float(l[1]), float(l[2]))] = np.array([float(l[3]),
                float(l[4]),
                float(l[5]),float(l[6])])

        # Line format of LAMMPS file -- > chunk x y z ncount vx vy vz
        with open(lammps_fname, "r") as lammps_file:
            lammps_lines = lammps_file.readlines()
        header = 3
        step_header = 1
        per_step = int(lammps_lines[3].split(" ")[1])
        begin = header + (step_header + per_step) * (s) + 1
        end = begin + per_step
        lammps_lines = lammps_lines[begin:end]
        lammps_lines = [l[:-1].split(" ") for l in lammps_lines]
        lammps_cells = {}
        for l in lammps_lines:
            l = filter(None, l)
            lammps_cells[(float(l[1]), float(l[2]), float(l[3]))] = np.array([float(l[5]),
                float(l[6]),
                float(l[7]),float(l[8])])
        for k in cfd_cells.keys():
            #print "CFD cell: " + str(cfd_cells[k])
            #print "LAMMPS cell: " + str(lammps_cells[k])
            try:
                if bin_dimension == 1:
                    lammps_cell = lammps_cells.values()[0]
                elif bin_dimension == 3:
                    lammps_cell = lammps_cells[k]
                diff_vel = abs(cfd_cells[k] - lammps_cell)
                if (np.any(diff_vel > tol)):
                    print "Step: %s" % s
                    print "Cell %s value differs in md : %s and cfd: %s" % (str(k), str(lammps_cell), str(cfd_cells[k]))
                    print "FAILURE"
                    if mode == "test":
                        assert False
                    else:
                        sys.exit()
            except KeyError:
                print "Step: %s" % s
                print "Cell not found :" + str(k)
                print "FAILURE"
                if mode == "test":
                    assert False
                else:
                    sys.exit()

    print "SUCCESS"
