#!/usr/bin/env python2
import sys, os
sys.path.append(os.path.abspath('..'))
from test_vels import compare_vels
       
bin_dimension = 1
if __name__ == "__main__":
    compare_vels(1e-6, steps=3, mode="debug", bin_dimension=bin_dimension)
