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

    "Socket" class for interfacing with CPL-Library.

Author(s)

    Eduardo Ramos

*/

#include <iostream>
#include "modify.h"
#include "domain.h"
#include "universe.h"
#include "input.h"
#include "comm.h"
#include "update.h"
#include "CPLSocketLAMMPS.h"


void CPLSocketLAMMPS::setTimingInfo() {
    initialStep = lmp->update->firststep;
    dt = lmp->update->dt;
}

void CPLSocketLAMMPS::setCartCommInfo() {
    int *npxyz = lmp->comm->procgrid;
    procGrid = std::vector<int>({npxyz[0], npxyz[1], npxyz[2]});
    myProcCoords = std::vector<int>({lmp->comm->myloc[0],
                                     lmp->comm->myloc[1],
                                     lmp->comm->myloc[2]});
}

void CPLSocketLAMMPS::setRealmDomainInfo() {
    double *global_domain = lmp->domain->prd;
    std::valarray<double> domain_orig({0.0 ,0.0, 0.0});
    std::valarray<double> domain_length({global_domain[0], 
                                         global_domain[1], 
                                         global_domain[2]});
    realmDomain = CPL::Domain(domain_orig, domain_length);
}


void CPLSocketLAMMPS::init() {
    CPLSocket::init();
	std::string cmd =  "variable CPLSTEPS equal " + std::to_string(nSteps);
    lmp->input->one(cmd.c_str());
}


void CPLSocketLAMMPS::configureBc(int mode) {
    double shift;
	if (mode == AVG_MODE_MIDPLANE) {
		shift = -(bcRegion.ly)/2.0;
	}
	else if (mode == AVG_MODE_BELOW) {
		shift = -bcRegion.ly;
	}
    // Create a new field with the corrected BC domain 
    std::valarray<double> domain_bounds = bcRegion.bounds;
    domain_bounds[2] += shift;
    domain_bounds[3] += shift;
    bcRegion = CPL::PortionField(CPL::Domain(domain_bounds), bcRegion.nCells,
                                 bcRegion.cellBounds);
    domain_bounds = bcPortionRegion.bounds;
    domain_bounds[2] += shift;
    domain_bounds[3] += shift;
    bcPortionRegion = CPL::PortionField(CPL::Domain(domain_bounds),
                                        bcPortionRegion.nCells, bcPortionRegion.cellBounds);
}
