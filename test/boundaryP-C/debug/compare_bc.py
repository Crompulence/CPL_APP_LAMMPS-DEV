#!/usr/bin/env python2
import sys, os
sys.path.append(os.path.abspath('..'))
from test_bc import compare_bc
       
bin_dimension = 1
if __name__ == "__main__":
    compare_bc(1e-5, steps=3, mode="debug", bin_dimension=bin_dimension)
