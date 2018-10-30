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

#include "fix_cpl_force.h"
#include "cpl/CPL_misclib.h"
#include "cpl/CPL_ndArray.h"

using namespace std::chrono;

FixCPLForce::FixCPLForce ( LAMMPS_NS::LAMMPS *lammps, int narg, char **arg) 
    : Fix (lammps, narg, arg)
{
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
                            | forceType_arg.compare("bndryavg") == 0)
                            break;
                        //Otherwise it is a sendtype argument and should be added
                        std::string forceType(*forcetype);
                        std::cout << "Lammps FixCPLForce forcetype: "  << forceType << " with args "
                                  <<  forceType_arg << std::endl;
                        forcetype_args.push_back(forcetype_arg);
                    }
                }
            }
        }
    }

    //int ifix = lammps->modify->find_fix("clumps");
    std::vector<double> fi{3};

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
    apply(1);
}


//NOTE -- Not actually called post force, for some reason
// this no longer works reliably in LAMMPS, instead call
// explicitly in CPLInit!
void FixCPLForce::apply(int nevery) {

    bool time = false;
    high_resolution_clock::time_point begin;
    high_resolution_clock::time_point end;

    if (time) begin = high_resolution_clock::now();

    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double *rmass = atom->rmass;
    double *radius = atom->radius;

    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    // This will get the total number of molecules including halos
    // However, acceleration won't be held for ghost cells so 
    // probably not possible to use this. Instead collect pre force
    // values and exchange these using CPL_swaphalos.
    int nall = atom->nlocal + atom->nghost;

    char* groupstr = "cplforcegroup";
    char* regionstr = "cplforceregion";

    int cplforcegroup = group->find(groupstr);
    int groupbit = group->bitmask[cplforcegroup];

    int rid = domain->find_region (regionstr);
    auto cplforceregion = domain->regions[rid];

    //Update CFD field buffer with latest recieved value
    fxyz->set_field(*cfdBuf);

    double fx=0., fy=0., fz=0.;
    double mi, radi, pot, xi[3], vi[3], ai[3];
    pot = 1.0; //Interaction Potential should be set here

    if (time) {
        end = high_resolution_clock::now();
        std::cout << " step " << update->ntimestep << " time allocation = " 
                 << duration_cast<microseconds>( end - begin ).count() << "e-6 s"   << std::endl;
        begin = high_resolution_clock::now();
    }

    //Only recalculate preforce everytime we recieve data
    if ((update->ntimestep%nevery == 0) | (fxyz->calc_preforce_everytime))
    {

        //Should we reset sums here?
        fxyz->resetsums();

        //Pre-force calculation, get quantities from discrete system needed to apply force
        if (fxyz->calc_preforce) {
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

    if (time) {
        end = high_resolution_clock::now();
        std::cout << " step " << update->ntimestep << " time pre force = "
                 << duration_cast<microseconds>( end - begin ).count() << "e-6 s"   << std::endl;
        begin = high_resolution_clock::now();
    }

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
            fi = fxyz->get_force(xi, vi, ai, mi, radi, pot);

            //Apply force and multiply by conversion factor if not SI or LJ units
            for (int n=0; n<3; n++){
                f[i][n] += fi[n]*force->ftm2v;
            }

//            std::cout.precision(17);
//            std::cout << "Force " << i << " " << xi[2] <<  " " << vi[2] << " " << ai[2] << " " <<
//                          mi << " " << fi[0] << " " << fi[1] << " " << fi[2] << " "  
//                          << f[i][0] << " " << f[i][1] << " " << f[i][2] << std::endl;

        }
    }

    if (time) {
        end = high_resolution_clock::now();
        std::cout <<  " step " << update->ntimestep << " time get force = " 
               << duration_cast<microseconds>( end - begin ).count() << "e-6 s"   << std::endl;
    }

}


void FixCPLForce::setupBuf (CPL::ndArray<double>& Buf, std::vector<int>& portion) {
    cfdBuf = &Buf;
	updateProcPortion(portion);
}
    

void FixCPLForce::updateProcPortion (std::vector<int>& portion) {

    procPortion.resize(6);
    for (int i = 0; i < 6; ++i) {
        procPortion[i] = portion[i];
    }
}
