#!/usr/bin/env python2
import sys, os
sys.path.append(os.path.abspath('..'))
from testutils import check_forces

check_forces(1e-6, mode="debug")
