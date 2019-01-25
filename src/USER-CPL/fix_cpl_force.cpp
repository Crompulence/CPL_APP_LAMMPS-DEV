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

    fix for creating forces to be applied as part of coupling with CPL-Library.

    This fix is turned on by fix_cpl_init but is also designed to potentiallt be used
    by itself for testing or other purposes. The usage looks like this:

    fix ID group-ID cpl/force region all forcetype X sendtype Y

    cpl/force returns a 3 per-atom values for the force which can be output using the dump custom 
    command, time-averaged with fix ave/atom command or reduced by the compute reduce command.
    Or the per-atom values can be referenced in an atom-style variable.
    
    An Example which prints two atom values:

    fix cplforcefix all cpl/force region all forcetype Drag peratom
    
    variable 	    cpl_out1 equal f_cplforcefix[1][3] 
    variable 	    cpl_out2 equal f_cplforcefix[4][1] 

    fix             print_stress all print 10 "${cpl_out1} ${cpl_out2}" file print_stress.txt screen yes
    dump            10 all custom 10000 dumpmyforce% id type x y z f_cplforcefix[*][*]

Author(s)

    David Trevelyan, Edward Smith, Eduardo Fernandez-Ramos

*/

#include<iostream>
#include<memory>
#include<fstream>
#include <cmath>
#include <chrono>

#include "atom.h"
#include "universe.h"
#include "error.h"
#include "group.h"
#include "domain.h"
#include "region.h"
#include "force.h"
#include "update.h"
#include "memory.h"

#include "fix_cpl_force.h"
#include "cpl/CPL_misclib.h"
#include "cpl/CPL_ndArray.h"
//#include "cpl/CPL_misclib.h"

using namespace std::chrono;


FixCPLForce::FixCPLForce ( LAMMPS_NS::LAMMPS *lammps, int narg, char **arg) 
    : Fix (lammps, narg, arg), fddata(NULL)
{

    calcperatom = false;
    class LAMMPS_NS::LAMMPS *lmp=lammps;
    for (int iarg=0; iarg<narg; iarg+=1){
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
                            | forceType_arg.compare("bndryavg") == 0) {
                            break;
                        } else if (forceType_arg.compare("calcperatom") == 0) {
                            calcperatom = true;
                            std::cout << " Output for fix cpl/force is on and can be accessed using f_cplforcefix[*][*] " << std::endl;
                        } else {

                            //Otherwise it is a forcetype argument and should be added
                            std::string forceType(*forcetype);
                            std::cout << "Lammps FixCPLForce forcetype: "  << forceType << " with args "
                                      <<  forceType_arg << std::endl;
                            forcetype_args.push_back(forcetype_arg);
                        }
                    }
                }
            }
        }
    }

    //int ifix = lammps->modify->find_fix("clumps");
    std::vector<double> fi{3};

    //Create peratom output
    if (calcperatom) {
        peratom_flag = 1;
        size_peratom_cols = numcols;
        peratom_freq = 1;

        grow_arrays(atom->nmax);
        //Set data array to zero
        int nlocal = atom->nlocal;
        for (int i = 0; i < nlocal; i++) {
        for (int j = 0; j < numcols; j++) {
            fddata[i][j] = 0; //i*3 + j;
        } }
    }
}

/* ---------------------------------------------------------------------- */

FixCPLForce::~FixCPLForce()
{
    if (calcperatom) {
        atom->delete_callback(id,0); //~ Unregister callbacks from the Atom class
        memory->destroy(fddata); //~ Destroy the local data array
    }
}


//NOTE: It is called from fixCPLInit initial_integrate() now.
int FixCPLForce::setmask() {
  int mask = 0;
  //mask |= LAMMPS_NS::FixConst::POST_FORCE;
  return mask;
}



/* ---------------------------------------------------------------------- */

void FixCPLForce::setup(int vflag)
{

    //Setup a map of default arguments for force types
    std::map <std::string, std::string> args_map{};
    for (int i=0; i<forcetype_args.size(); i=i+2) {
        args_map.insert( { *forcetype_args[i], *forcetype_args[i+1]});
    }

    //This is a factory
    std::string fxyzType(*forcetype);
    if (fxyzType.compare("Flekkoy") == 0) {
        fxyz = std::make_unique<CPLForceFlekkoy>(cfdBuf->shape(0), 
                                                 cfdBuf->shape(1), 
                                                 cfdBuf->shape(2), 
                                                 cfdBuf->shape(3));
    } else if (fxyzType.compare("test") == 0) {
        fxyz = std::make_unique<CPLForceTest>(cfdBuf->shape(0), 
                                              cfdBuf->shape(1), 
                                              cfdBuf->shape(2), 
                                              cfdBuf->shape(3));
    } else if (fxyzType.compare("Velocity") == 0) {
        fxyz = std::make_unique<CPLForceVelocity>(cfdBuf->shape(0), 
                                                  cfdBuf->shape(1), 
                                                  cfdBuf->shape(2), 
                                                  cfdBuf->shape(3));
    } else if (fxyzType.compare("Drag") == 0) {
        fxyz = std::make_unique<CPLForceDrag>(cfdBuf->shape(0), 
                                              cfdBuf->shape(1), 
                                              cfdBuf->shape(2), 
                                              cfdBuf->shape(3), 
                                              args_map);
        //fxyz->calc_preforce = 1;
    } else if (fxyzType.compare("Di_Felice") == 0) {
        fxyz = std::make_unique<CPLForceGranular>(cfdBuf->shape(0), 
                                                  cfdBuf->shape(1), 
                                                  cfdBuf->shape(2), 
                                                  cfdBuf->shape(3), 
                                                  args_map);  

    } else if (fxyzType.compare("Ergun") == 0) {
        fxyz = std::make_unique<CPLForceErgun>(cfdBuf->shape(0), 
                                               cfdBuf->shape(1), 
                                               cfdBuf->shape(2), 
                                               cfdBuf->shape(3), 
                                               args_map); 

    } else if (fxyzType.compare("BVK") == 0) {
        fxyz = std::make_unique<CPLForceBVK>(cfdBuf->shape(0), 
                                             cfdBuf->shape(1), 
                                             cfdBuf->shape(2), 
                                             cfdBuf->shape(3), 
                                             args_map); 
    } else {
        std::string cmd("CPLForce type ");
        cmd += fxyzType + " not defined";
        throw std::runtime_error(cmd);
    }

    //Set CPLForce min/max to local processor limits using values from CPL library 
	double min[3]; double max[3];
    std::vector<int> cnstFPortion(6);
    std::vector<int> cnstFRegion(6);
    CPL::get_cnst_limits(cnstFRegion.data());
    CPL::my_proc_portion (cnstFRegion.data(), cnstFPortion.data());
	//MIN
	CPL::map_cell2coord(cnstFPortion[0], 
                        cnstFPortion[2], 
                        cnstFPortion[4], min);
	//MAX
	CPL::map_cell2coord(cnstFPortion[1], 
                        cnstFPortion[3], 
                        cnstFPortion[5], max);

    double dx = CPL::get<double> ("dx");     
    double dy = CPL::get<double> ("dy");     
    double dz = CPL::get<double> ("dz");     

	max[0] += dx;
	max[1] += dy;
	max[2] += dz;
    fxyz->set_minmax(min, max);

    //Call apply for first step
    //apply(1, 1, 1);
    irepeat = 0;
}



void FixCPLForce::pre_force(int Nfreq, int Nrepeat, int Nevery){

    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double *rmass = atom->rmass;
    double *radius = atom->radius;

    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double mi, radi, pot, xi[3], vi[3], ai[3];
    pot = 1.0; //Interaction Potential should be set here

    char* groupstr = "cplforcegroup";
    int cplforcegroup = group->find(groupstr);
    int groupbit = group->bitmask[cplforcegroup];

    char* regionstr = "cplforceregion";
    int rid = domain->find_region (regionstr);
    auto cplforceregion = domain->regions[rid];

    //Update CFD field buffer with latest recieved value
    fxyz->set_field(*cfdBuf);

    //Only recalculate preforce everytime we recieve data
    // or Nevery as this accumulates data for send as required
    if ((update->ntimestep%Nevery == 0) | (fxyz->calc_preforce_everytime))
    {

         //If Nrepeat, then reset sums
         if (irepeat == Nrepeat){
#if DEBUG
            std::cout <<  "Resetting sums " <<  irepeat << " " << Nrepeat << std::endl;
#endif
            reset_sums();
         } else {
            irepeat++;
#if DEBUG
            std::cout <<  "irepeat = " <<  irepeat << " of " << Nrepeat << std::endl;
#endif
         }

        //Pre-force calculation, get quantities from discrete system needed to apply force
        if (fxyz->calc_preforce) {

            //Increment pre-force counter
            fxyz->Npre_force++;

            //All local particles
        	for (int i = 0; i < nlocal; ++i)
        	{
           		if (mask[i] & groupbit)
            	{
		            //Get local molecule data
		            if (atom->rmass_flag) {mi = rmass[i];}
                    else {mi = 1;}
		            if (atom->radius_flag) {radi = radius[i];}
                    else {radi = 1;}
		            for (int n=0; n<3; n++){
		                xi[n]=x[i][n]; 
		                vi[n]=v[i][n]; 
		                ai[n]=f[i][n];
		            }

		            // Sum all the weights for each cell.
		            fxyz->pre_force(xi, vi, ai, mi, radi, pot);

            	}
            }
        }
    }

}


void FixCPLForce::apply_force(int Nfreq, int Nrepeat, int Nevery){

    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double *rmass = atom->rmass;
    double *radius = atom->radius;

    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double mi, radi, pot, xi[3], vi[3], ai[3];
    pot = 1.0; //Interaction Potential should be set here

    char* groupstr = "cplforcegroup";
    int cplforcegroup = group->find(groupstr);
    int groupbit = group->bitmask[cplforcegroup];

    char* regionstr = "cplforceregion";
    int rid = domain->find_region (regionstr);
    auto cplforceregion = domain->regions[rid];

    //std::cout << "pre set field " << fxyz->Nforce << " " << cfdBuf->shape(0) << " " << 
    //          cfdBuf->shape(1) << " " << cfdBuf->shape(2) << " " << cfdBuf->shape(3) << std::endl;

    //Update CFD field buffer with latest recieved value
    fxyz->set_field(*cfdBuf);

    //Increment force counter
    fxyz->Nforce++;

    // Calculate force and apply
    for (int i = 0; i < nlocal; ++i)
    {
        if (mask[i] & groupbit)
        {

            //Get local molecule data
            if (atom->rmass_flag) {mi = rmass[i];}
            else {mi = 1;}
            if (atom->radius_flag) {radi = radius[i];}
            else {radi = 1;}
            for (int n=0; n<3; n++){
                xi[n]=x[i][n]; 
                vi[n]=v[i][n]; 
                ai[n]=f[i][n];
            }

            //Get force from object
            //std::vector<double> fi = {0., 0., 0.};
            auto fi = fxyz->get_force(xi, vi, ai, mi, radi, pot);

            //Apply force and multiply by conversion factor if not SI or LJ units
            for (int n=0; n<3; n++){
                f[i][n] += fi[n]*force->ftm2v;
                if (calcperatom) fddata[i][n] = fi[n]*force->ftm2v;
            }

            //std::cout.precision(17);
//            std::cout << "Force " <<  update->ntimestep << " " << i << " " << mi << " " 
//                      << xi[1] <<  " " << vi[0] << " " << ai[0] << " " << fi[0] << " " 
//                          << f[i][0] << " " << f[i][1] << " " << f[i][2] << std::endl;

        }
    }

}


void FixCPLForce::post_constraint_force(int Nfreq, int Nrepeat, int Nevery){

    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double *rmass = atom->rmass;
    double *radius = atom->radius;

    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double mi, radi, pot, xi[3], vi[3], ai[3];
    pot = 1.0; //Interaction Potential should be set here

    char* groupstr = "cplforcegroup";
    int cplforcegroup = group->find(groupstr);
    int groupbit = group->bitmask[cplforcegroup];

    char* regionstr = "cplforceregion";
    int rid = domain->find_region (regionstr);
    auto cplforceregion = domain->regions[rid];

    // Only recalculate post force everytime we recieve data
    // or Nevery as this accumulates data for send as required
    if ((update->ntimestep%Nevery == 0)) //| (fxyz->calc_postforce_everytime))
    {

        //We cannot reset sums here as we would lose Fsum

        //Post-force calculation, get quantities from discrete system needed to apply force
        if (fxyz->calc_postforce) {

            //Increment post force counter
            fxyz->Npost_force++;

            //All local atoms
        	for (int i = 0; i < nlocal; ++i)
        	{
           		if (mask[i] & groupbit)
            	{
	                //Get local molecule data
	                if (atom->rmass_flag) {mi = rmass[i];}
                    else {mi = 1;}
	                if (atom->radius_flag) {radi = radius[i];}
                    else {radi = 1;}
	                for (int n=0; n<3; n++){
	                    xi[n]=x[i][n]; 
	                    vi[n]=v[i][n]; 
	                    ai[n]=f[i][n];
	                }

	                // Sum all the weights for each cell.
	                fxyz->post_force(xi, vi, ai, mi, radi, pot);

            	}
            }
        }

    }

}

//NOTE -- Not actually called post force, for some reason
// this no longer works reliably in LAMMPS, instead call
// explicitly in CPLInit!
void FixCPLForce::apply(int Nfreq, int Nrepeat, int Nevery) {

    bool time = false;
    high_resolution_clock::time_point begin;
    high_resolution_clock::time_point end;

    if (time) begin = high_resolution_clock::now();

//    char* groupstr = "cplforcegroup";
//    char* regionstr = "cplforceregion";

//    int cplforcegroup = group->find(groupstr);
//    int groupbit = group->bitmask[cplforcegroup];

//    int rid = domain->find_region (regionstr);
//    auto cplforceregion = domain->regions[rid];

    if (time) {
        end = high_resolution_clock::now();
        std::cout << " step " << update->ntimestep << " time allocation = " 
                 << duration_cast<microseconds>( end - begin ).count() << "e-6 s"   << std::endl;
        begin = high_resolution_clock::now();
    }

    // Do calculations required before applying force
    pre_force(Nfreq, Nrepeat, Nevery);

    if (time) {
        end = high_resolution_clock::now();
        std::cout << " step " << update->ntimestep << " time pre force = "
                 << duration_cast<microseconds>( end - begin ).count() << "e-6 s"   << std::endl;
        begin = high_resolution_clock::now();
    }

    //Apply force
    apply_force(Nfreq, Nrepeat, Nevery);

    if (time) {
        end = high_resolution_clock::now();
        std::cout <<  " step " << update->ntimestep << " time get force = " 
               << duration_cast<microseconds>( end - begin ).count() << "e-6 s"   << std::endl;
    }

}


void FixCPLForce::setupBuf(CPL::ndArray<double>& Buf, std::vector<int>& portion) {
    cfdBuf = &Buf;
	updateProcPortion(portion);
}
    

void FixCPLForce::updateProcPortion (std::vector<int>& portion) {

    procPortion.resize(6);
    for (int i = 0; i < 6; ++i) {
        procPortion[i] = portion[i];
    }
}


void FixCPLForce::reset_sums(){
#if DEBUG
     std::cout <<  "Resetting sums " <<  irepeat << std::endl;
#endif
     fxyz->resetsums();
     irepeat = 0;
}



/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixCPLForce::grow_arrays(int nmax)
{
  memory->grow(fddata,nmax,numcols,"cplforce:fddata");
  array_atom = fddata;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixCPLForce::copy_arrays(int i, int j, int delflag)
{
  for (int q = 0; q < numcols; q++)
    fddata[j][q] = fddata[i][q];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixCPLForce::pack_exchange(int i, double *buf)
{
  for (int q = 0; q < numcols; q++)
    buf[q] = fddata[i][q];

  return numcols;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixCPLForce::unpack_exchange(int nlocal, double *buf)
{
  for (int q = 0; q < numcols; q++)
    fddata[nlocal][q] = buf[q];

  return numcols;
}

// This will get the total number of molecules including halos
// However, acceleration won't be held for ghost cells so 
// probably not possible to use this. Instead collect pre force
// values and exchange these using CPL_swaphalos.
//int nall = atom->nlocal + atom->nghost;
