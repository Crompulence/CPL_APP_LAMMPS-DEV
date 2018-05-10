#!/usr/bin/env python
import os
import sys
import glob
from shutil import copyfile
from datetime import datetime

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

DIR = os.path.dirname(os.path.realpath(__file__))
with cd(DIR):
    os.remove("./mpmd.patch")
    patches = glob.glob("./*.patch")

    lammpssrcdir = sys.argv[1]
    with open(lammpssrcdir + "/version.h") as f:
        #Get datetime of lammps
        s = f.read()
        l = s.replace('"','').split()
        lammps_version_date = datetime.strptime(l[-2] +" " + l [-3] + " " + l[-1], '%b %d %Y')

        #Get datetime of patches
        patch_to_apply = datetime.strptime("Jan 30 2100", '%b %d %Y') 
        for p in patches:
            pd = p.replace("./","").replace("mpmd_","").replace(".patch","")
            patchdate = datetime.strptime(pd[-8:-4] + " " + str(lammps_version_date.day) +" " + pd[-4:], '%b %d %Y')
            #print(patchdate, lammps_version_date, patchdate > lammps_version_date)
            if patchdate >= lammps_version_date:
                #Use earliest possible patch
                if patch_to_apply > patchdate:
                    patch_to_apply = patchdate
                    patchfile = p
                else:
                    pass


        print(s, " applying patch ", patchfile)

        #Take the latest one and use this
        copyfile(patchfile, "./mpmd.patch")




