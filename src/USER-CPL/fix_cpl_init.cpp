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

           Copyright (C) 2012-2018 Edward Smith & David Trevelyan

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



#include<iostream>
#include <string.h>
#include <stdlib.h>

#include "update.h"
#include "error.h"

#include "fix_cpl_init.h"

fixCPLInit::fixCPLInit(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    		: Fix (lammps, narg, arg) 
{
    class LAMMPS_NS::LAMMPS *lmp=lammps;
    cplsocket.initComms();
    cplsocket.initMD(lmp);
    nevery = cplsocket.timestep_ratio;

    forcetype = std::make_shared<std::string>("Undefined");
    sendtype = std::make_shared<std::string>("velocity");
    std::vector<std::shared_ptr<std::string>> sendtype_list;
    bndryavg = std::make_shared<std::string>("above");     //default to above if not specified


    for (int iarg=0; iarg<narg; iarg+=1){
        std::cout << "Lammps cpl/init input arg "  << iarg << " is " << arg[iarg] << std::endl;
        std::string arguments(arg[iarg]);
        if (arguments == "forcetype"){
            if (iarg+1<narg) {
                forcetype = std::make_shared<std::string>(arg[iarg+1]);
                if (iarg+2<narg) {
                    for (int jarg=iarg+2; jarg<narg; jarg+=1){
                        std::shared_ptr<std::string> forcetype_arg;
                        forcetype_arg = std::make_shared<std::string>(arg[jarg]);
                        //Check if we have read another argument type
                        std::string forceType_arg(*forcetype_arg);
                        if (  forceType_arg.compare("sendtype") == 0 
                            | forceType_arg.compare("bndryavg") == 0)
                            break;
                        //Otherwise it is a sendtype argument and should be added
                        std::string forceType(*forcetype);
                        std::cout << "Lammps forcetype: "  << forceType << " with args "
                                  <<  forceType_arg << std::endl;
                        forcetype_args.push_back(forcetype_arg);
                    }
                }
            }
        }

        if (arguments == "sendtype")
            if (iarg+1<narg) {
                for (int jarg=iarg+1; jarg<narg; jarg+=1){
                    sendtype = std::make_shared<std::string>(arg[jarg]);
                    //Check if we have read another argument type
                    std::string sendType(*sendtype);
                    if (  sendType.compare("forcetype") == 0 
                        | sendType.compare("bndryavg") == 0)
                        break;
                    //Otherwise it is a sendtype argument and should be added
                    sendtype_list.push_back(sendtype);
                }
            }

        if (arguments == "bndryavg"){
            if (iarg+1<narg)
                bndryavg = std::make_shared<std::string>(arg[iarg+1]);
        }

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
    sendbitflag = 0;
    for ( auto &sendtype : sendtype_list ) {
        std::string sendType(*sendtype);
        //Pick 'n' mix send types
        if (sendType.compare("VEL") == 0){
            sendbitflag = sendbitflag | cplsocket.VEL;
        } else if (sendType.compare("NBIN") == 0) {
            sendbitflag = sendbitflag | cplsocket.NBIN;
        } else if (sendType.compare("STRESS") == 0) {
            sendbitflag = sendbitflag | cplsocket.STRESS;
        } else if (sendType.compare("FORCE") == 0) {
            sendbitflag = sendbitflag | cplsocket.FORCE;
        } else if (sendType.compare("FORCECOEFF") == 0) {
            sendbitflag = sendbitflag | cplsocket.FORCECOEFF;
        } else if (sendType.compare("VOIDRATIO") == 0) {
            sendbitflag = sendbitflag | cplsocket.VOIDRATIO;
        //Predefined sendtypes
        } else if (sendType.compare("velocity") == 0) {
            sendbitflag = sendbitflag | cplsocket.VEL | cplsocket.NBIN;
        } else if (sendType.compare("gran") == 0) {
            sendbitflag = sendbitflag | cplsocket.FORCE | cplsocket.VOIDRATIO;
        } else if (sendType.compare("granfull") == 0) {
            sendbitflag = sendbitflag | cplsocket.VEL | cplsocket.FORCE |
                          cplsocket.FORCECOEFF | cplsocket.VOIDRATIO;
        } else { 
            std::cout << "Lammps sendtype: "  << sendType << " not recognised"
                      << " bitflag so far: " << sendbitflag << std::endl;
            lammps->error->all(FLERR,"Lammps sendtype Error");
        }

    }

//    if (((sendbitflag & cplsocket.FORCECOEFF) == cplsocket.FORCECOEFF)
//         & (forceType.compare("Drag") != 0)) {
//        lammps->error->all(FLERR,"Drag Forcetype (or its derivatives) required for sendtype granfull");
//    }

}

void fixCPLInit::init(){

}

void fixCPLInit::setas_last_fix() {
   int ifix = lmp->modify->find_fix("cplfix");
   LAMMPS_NS::Fix* fix_aux = lmp->modify->fix[ifix];
   int nfix = lmp->modify->nfix; 
   lmp->modify->fix[ifix] = lmp->modify->fix[nfix-1];
   lmp->modify->fix[nfix-1] = fix_aux;
   lmp->modify->fmask[ifix] = lmp->modify->fix[ifix]->setmask();
   lmp->modify->fmask[nfix-1] = lmp->modify->fix[nfix-1]->setmask();
}


int fixCPLInit::setmask() {
  int mask = 0;
  mask |= LAMMPS_NS::FixConst::POST_FORCE;
  return mask;
}

void fixCPLInit::post_constructor() {
	
    //Setup what to send and how to apply forces
    cplsocket.setupFixMDtoCFD(lmp, sendbitflag);

    //Note that constraint fix is setup through lammps input system
    cplsocket.setupFixCFDtoMD(lmp, forcetype, forcetype_args);

}


void fixCPLInit::setup(int vflag)
{
  	post_force(vflag);
}



void fixCPLInit::post_force(int vflag)
{

    //std::cout << "fixCPLInit " << update->ntimestep << " " << update->ntimestep%nevery <<" " <<  update->dt << " " << std::endl;   

    // Recieve and unpack from CFD
    if (update->ntimestep%nevery == 0){
        cplsocket.receive();
    }
    cplsocket.cplfix->apply(nevery);

    //Pack and send to CFD
    if (update->ntimestep%nevery == 0){
        cplsocket.pack(lmp, sendbitflag);
        cplsocket.send();
    }


}

fixCPLInit::~fixCPLInit() {
	cplsocket.finalizeComms();
}


