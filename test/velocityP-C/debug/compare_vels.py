#!/usr/bin/env python2
import sys, os
sys.path.append(os.path.abspath('..'))
from test_vels import check_vels
       

if __name__ == "__main__":
    check_vels(1e-6, steps=3, mode="debug")
