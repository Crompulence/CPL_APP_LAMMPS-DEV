
# Bouncing Ball Test

This directory is the test for a single particle driven by an external force, 
which is recieved from a coupled mock using CPL library,
towards a wall, before bouncing on the wall.

## What does it do?

It runs the case of a single particle by running
a coupled case with lammps and a mock CFD script
before checking the solution against the constant 
accelerations equations for the falling bit before/after it 
hits the wall and using regression for the wall interaction.

![alt tag](https://raw.githubusercontent.com/Crompulence/CPL_APP_LAMMPS-DEV/blob/master/test/constant_force/Wall_force_LAMMPS.png)

This is run for a range of processor topologies by pytest and checked each time
against the constant acceleration equations.

```python
@pytest.mark.parametrize("mdprocs", [1, 2, 4, 8])
def test_gravity(build_run, mdprocs):

    clean_dir()
    run = run_case(mdprocs)

    #Check vs analystical solution for gravity
    import bouncing
    with cd(TEST_DIR):
        error = bouncing.check_bouncing_error_vs_gravity()
        for e in error[0,1,:]:
            assert np.abs(e) < 1e-11
```

 Which is linked to the CFD mock script

```python
import numpy as np
from mpi4py import MPI
from cplpy import CPL

g = 9.81
mi = -5.9490638385009208e-08

#initialise MPI and CPL
comm = MPI.COMM_WORLD
CPL = CPL()
MD_COMM = CPL.init(CPL.CFD_REALM)

#Setup coupled simulation
CPL.setup_cfd(MD_COMM.Create_cart([1, 1, 1]), [1.5E-003, 1.5E-003, 2.5E-003],
             [0.0, 0.0, 0.0], [8, 8, 8])
recv_array, send_array = CPL.get_arrays(recv_size=4, send_size=9)

#Run 100 steps
for time in range(101):
    send_array[2,:,:,:] = mi*g
    CPL.send(send_array)
    recv_array, ierr = CPL.recv(recv_array)
    print("CFD time = ", time)

comm.Barrier()
CPL.finalize()
MPI.Finalize()
```

which can be 

## What is tested?

This tests the information exchange between the codes
by CPL library, application of thesimplest possible force to the particle,
updating of position based on the force including the moving of
 a particle over multiple processors in LAMMPS.

## If it breaks, what do you check

If all other tests in CPL library [https://travis-ci.org/Crompulence/cpl-library] are passing up to this point,
it's possible that CPL_force has been changed to break this. 
Alternativly
