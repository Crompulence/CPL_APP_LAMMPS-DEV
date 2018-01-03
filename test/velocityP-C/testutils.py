import numpy as np
import re
import sys

def name_with_step(fname, step):
    split_fname = re.split('(\W)', fname)
    dot_pos = split_fname.index(".")
    split_fname.insert(dot_pos , str(step))
    return "".join(split_fname)

def check_vels(tol, lammps_fname="lammps_vels.dat", cfd_fname="cfd_vels.dat", steps=1, mode="test"):

    # Line format of CFD script file -- > x y z vx vy vz
    for s in xrange(steps):
        with open(name_with_step(cfd_fname, s), "r") as cfd_file:
            cfd_lines = cfd_file.readlines()
        cfd_lines = [l[:-1].split(" ") for l in cfd_lines]
        cfd_cells = {}
        for l in cfd_lines:
            cfd_cells[(float(l[0]), float(l[1]), float(l[2]))] = np.array([float(l[3]),
                float(l[4]),
                float(l[5])])

        # Line format of LAMMPS file -- > chunk x y z ncount vx vy vz
        with open(lammps_fname, "r") as lammps_file:
            lammps_lines = lammps_file.readlines()
        header = 3
        step_header = 1
        per_step = int(lammps_lines[3].split(" ")[1])
        #NOTE: Jump the first time-step. Initial state.
        begin = header + (step_header + per_step) * (s+1) + 1
        end = begin + per_step
        lammps_lines = lammps_lines[begin:end]
        lammps_lines = [l[:-1].split(" ") for l in lammps_lines]
        lammps_cells = {}
        for l in lammps_lines:
            l = filter(None, l)
            lammps_cells[(float(l[1]), float(l[2]), float(l[3]))] = np.array([float(l[5]),
                float(l[6]),
                float(l[7])])
        for k in cfd_cells.keys():
            #print "CFD cell: " + str(cfd_cells[k])
            #print "LAMMPS cell: " + str(lammps_cells[k])
            try:
                diff_vel = abs(cfd_cells[k] - lammps_cells[k])
                if (np.any(diff_vel > tol)):
                    print "Step: %s" % s
                    print "Cell %s value differs in md : %s and cfd: %s" % (str(k), str(lammps_cells[k]), str(cfd_cells[k]))
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
 
