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



fixCPLInit::fixCPLInit(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    		: Fix (lammps, narg, arg) {

	if (narg < 8) error->all(FLERR,"Illegal fix cplinit command");
	int bndry_avg_mode;
	if (strcmp(arg[3],"below") == 0)
		bndry_avg_mode = AVG_MODE_BELOW;
	else if (strcmp(arg[3],"above") == 0)
		bndry_avg_mode = AVG_MODE_ABOVE;
	else if (strcmp(arg[3],"midplane") == 0)
		bndry_avg_mode = AVG_MODE_MIDPLANE;
	else
		error->all(FLERR,"Illegal fix cplinit command - averaging method should \
				   be 'below', 'above' or 'midplane' ");
    int units;
	if (strcmp(arg[5],"real") == 0)
        cplsocket.units = REAL_UNITS;
    else if (strcmp(arg[5],"lj") == 0)
        cplsocket.units = LJ_UNITS;
	else
		error->all(FLERR,"Illegal fix cplinit command - units should be 'real' or 'lj'");
   	class LAMMPS_NS::LAMMPS *lmp=lammps;
   	cplsocket.initMD(lammps);
	cplsocket.setBndryAvgMode(bndry_avg_mode);
   	nevery = cplsocket.timestep_ratio;
   	std::cout << "LAMMPS NEVERY: " << nevery <<std::endl;
}

int fixCPLInit::setmask() {
  int mask = 0;
  mask |= LAMMPS_NS::FixConst::END_OF_STEP;
  mask |= LAMMPS_NS::FixConst::POST_FORCE;
  return mask;
}


void fixCPLInit::init() {
	
   cplsocket.setupFixCFDtoMD(lmp);
   cplsocket.setupFixMDtoCFD(lmp);
}


void fixCPLInit::setup(int vflag) {
  	end_of_step();
}


void fixCPLInit::post_force(int vflag) {
	cplsocket.updateStress();
}

void fixCPLInit::end_of_step() {
	static int c = 0;
    // Communications
	//std::cout << "SEND: " <<  c << std::endl;
	c++;
    cplsocket.recvStress();
	cplsocket.cfdbcfix->end_of_step();
	cplsocket.packVelocity(lmp);
    cplsocket.sendVelocity();
}

fixCPLInit::~fixCPLInit() {
	cplsocket.finalizeComms();
}

