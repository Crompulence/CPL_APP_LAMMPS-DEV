/*

    ________/\\\\\\\\\__/\\\\\\\\\\\\\____/\\\_____________
     _____/\\\////////__\/\\\/////////\\\_\/\\\_____________
      ___/\\\/___________\/\\\_______\/\\\_\/\\\_____________
       __/\\\_____________\/\\\\\\\\\\\\\/__\/\\\_____________
        _\/\\\_____________\/\\\/////////____\/\\\_____________
         _\//\\\____________\/\\\_____________\/\\\_____________
          __\///\\\__________\/\\\_____________\/\\\_____________
           ____\////\\\\\\\\\_\/\\\_____________\/\\\\\\\\\\\\\\\_
            _______\/////////__\///______________\///////////////__


                         C P L  -  L I B R A R Y

           Copyright (C) 2012-2015 Edward Smith & David Trevelyan

License

    This file is part of CPL-Library.

    CPL-Library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CPL-Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CPL-Library.  If not, see <http://www.gnu.org/licenses/>.

Description

    "Initialiser fix" for coupled simulation with CPL-Library.

Author(s)

    Edward Smith, Eduardo Ramos Fernandez

*/



#include "fix_cpl_init.h"
#include<iostream>
#include <string.h>
#include "error.h"
#include <stdlib.h>

#include<iostream>
#include "fix_cpl_init.h"
#include "update.h"

fixCPLInit::fixCPLInit(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    		: Fix (lammps, narg, arg) 
{
    class LAMMPS_NS::LAMMPS *lmp=lammps;
    cplsocket.initComms();
    cplsocket.initMD(lmp);
    nevery = cplsocket.timestep_ratio;

    forcetype = std::make_shared<std::string>("Undefined");
    sendtype = std::make_shared<std::string>("velocity");
    //Change from undefined to default of above if not specified
    bndryavg = std::make_shared<std::string>("above");
    for (int iarg=0; iarg<narg; iarg+=1){
        std::cout << iarg << " " << arg[iarg] << std::endl;
        std::string arguments(arg[iarg]);
        if (arguments == "forcetype")
            if (iarg+1<narg)
                forcetype = std::make_shared<std::string>(arg[iarg+1]);

        if (arguments == "sendtype")
            if (iarg+1<narg)
                sendtype = std::make_shared<std::string>(arg[iarg+1]);

        if (arguments == "bndryavg")
            if (iarg+1<narg)
                bndryavg = std::make_shared<std::string>(arg[iarg+1]);

    }
    //Raise error if forcetype is not specified
    std::string forceType(*forcetype);
    if (forceType.compare("Undefined") == 0){
        lammps->error->all(FLERR,"Must specify forcetype on cpl/init line in LAMMPS input file");
    }

    //Raise error if bndry_avg is not specified, otherwise set this
    std::string BndryAvg(*bndryavg);
    //if (BndryAvg.compare("Undefined") == 0){
    //    lammps->error->all(FLERR,"Illegal fix cplinit command - bndryavg should be 'below', 'above' or 'midplane'");
    //} else if (BndryAvg.compare("below") == 0) {
    if (BndryAvg.compare("below") == 0) {
        cplsocket.setBndryAvgMode(AVG_MODE_BELOW);
    } else if (BndryAvg.compare("above") == 0) {
        cplsocket.setBndryAvgMode(AVG_MODE_ABOVE);
    } else if (BndryAvg.compare("midplane") == 0) {
    	cplsocket.setBndryAvgMode(AVG_MODE_MIDPLANE);
    }

    //Create appropriate bitflag to determine what is sent
    std::string sendType(*sendtype);
    if (sendType.compare("velocity") == 0){
        sendbitflag = cplsocket.VEL;
    } else if (sendType.compare("gran") == 0) {
        //cplsocket.packGran(lmp);
        sendbitflag = cplsocket.FORCE | cplsocket.VOIDRATIO;
    } else if (sendType.compare("granfull") == 0) {
        sendbitflag = cplsocket.VEL | cplsocket.FORCE |
                      cplsocket.FORCECOEFF | cplsocket.VOIDRATIO;
    }

}

int fixCPLInit::setmask() {
  int mask = 0;
  mask |= LAMMPS_NS::FixConst::POST_FORCE;
  return mask;
}


void fixCPLInit::init()
{
	
    //Setup what to send and how to apply forces
    cplsocket.setupFixMDtoCFD(lmp, sendbitflag);
    cplsocket.setupFixCFDtoMD(lmp, forcetype);

}


void fixCPLInit::setup(int vflag)
{
  	post_force(vflag);
}



void fixCPLInit::post_force(int vflag)
{

    //std::cout << "fixCPLInit " << update->ntimestep << " " << update->ntimestep%nevery << std::endl;   

    // Recieve and unpack from CFD
    if (update->ntimestep%nevery == 0){
        cplsocket.receive();
    }
    cplsocket.cplfix->apply();

    //Pack and send to CFD
    if (update->ntimestep%nevery == 0){
        cplsocket.pack(lmp, sendbitflag);
        cplsocket.send();
    }


}

fixCPLInit::~fixCPLInit() {
	cplsocket.finalizeComms();
}

