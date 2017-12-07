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

    "CPLSocketLAMMPS" class for interfacing with CPL-Library.

Author(s)

    Eduardo Ramos

*/
#ifndef CPL_SOCKET_LAMMPS_H_INCLUDED
#define CPL_SOCKET_LAMMPS_H_INCLUDED

#include <vector>
#include <valarray>
#include "mpi.h"
#include "lammps.h"
#include "fix_cpl_force.h"
#include "cpl/cpl.h"
#include "cpl/CPL_ndArray.h"
#include "cpl/CPL_field.h"
#include "cpl/CPLSocket.h"


const int AVG_MODE_ABOVE = 0;
const int AVG_MODE_BELOW = 1;
const int AVG_MODE_MIDPLANE = 2;

class CPLSocketLAMMPS : public CPLSocket
{

public:
    
    // Construct from no arguments
    CPLSocketLAMMPS() : CPLSocket(CPL::md_realm){}; 
    virtual ~CPLSocketLAMMPS(){};
    void init();

    // Data preparation and communication 
    void configureBc(int mode);
    void setTimingInfo();
    void setCartCommInfo();
    void setRealmDomainInfo();
    void setLammps(LAMMPS_NS::LAMMPS* lammps) {lmp = lammps;}

    //Specific attributes
    FixCPLForce* cplfix;
    LAMMPS_NS::LAMMPS* lmp;

};

#endif // CPL_SOCKET_LAMMPS_H_INCLUDED
