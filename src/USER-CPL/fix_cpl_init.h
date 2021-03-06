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

    Edward Smith, Eduardo Fernandez-Ramos, David Trevelyan

*/
#ifdef FIX_CLASS

FixStyle(cpl/init, fixCPLInit)

#else

#ifndef LMP_FIX_CPL_INIT_H
#define LMP_FIX_CPL_INIT_H

#include <memory>

#include "fix.h"
#include "modify.h"

#include "cpl/cpl.h"
//#include "cpl/TransmittingField.h"

#include "CPLSocketLAMMPS.h"


class fixCPLInit : public LAMMPS_NS::Fix {

public:

    fixCPLInit(class LAMMPS_NS::LAMMPS *lammps, int narg, char **arg);
	~fixCPLInit();
    int setmask();
    void init (); 
    void setup (int vflag); 
	void post_integrate();
	void post_force(int vflag);
    void post_constructor();
    void setas_last_fix();

    // Pool of boundary conditions for CFD
    //CPL::OutgoingFieldPool bcPool;
    //CPL::IncomingFieldPool cnstPool;

    CPLSocketLAMMPS cplsocket;
    std::shared_ptr<std::string> forcetype;
    std::shared_ptr<std::string> sendtype;
    std::shared_ptr<std::string> bndryavg;
    int sendbitflag;
    std::vector<std::shared_ptr<std::string>> forcetype_args;
    //Send averaging frequencies
    int Nfreq, Nrepeat, Nevery;

private:
    LAMMPS_NS::LAMMPS* lmp;

};

#endif
#endif
