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

    "Socket" class for interfacing with CPL-Library.

Author(s)

    David Trevelyan, Edward Smith, Eduardo Fernandez-Ramos

*/
#include<iostream>
#include <iomanip>

#include "update.h"
#include "modify.h"
#include "fix_ave_chunk.h"
#include "domain.h"
#include "universe.h"
#include "input.h"
#include "comm.h"
#include "error.h"

#include "CPLSocketLAMMPS.h"
#include "cpl/CPL_cartCreate.h"


void CPLSocketLAMMPS::initComms() {

    // Split MPI_COMM_WORLD into realm communicators
    CPL::init(CPL::md_realm, realmComm);
    MPI_Comm_rank(realmComm, &rankRealm);

};

void CPLSocketLAMMPS::finalizeComms() {
    CPL::finalize();
};

void CPLSocketLAMMPS::initMD(LAMMPS_NS::LAMMPS *lammps) {

    // Store my own coordinates for later
    myCoords[0] = lammps->comm->myloc[0];
    myCoords[1] = lammps->comm->myloc[1];
    myCoords[2] = lammps->comm->myloc[2];

    // Parameters for coupler_md_init
    int initialstep = lammps->update->firststep;
    double dt = lammps->update->dt;
    int *npxyz_md = lammps->comm->procgrid;
    double *globaldomain = lammps->domain->prd;
    double dummydensity = -666.0;

    // Set up new cartesian communicator with same coordinates as lammps
    // interal cartesian communicator (based on mycoords)
    MPI_Comm icomm_grid;
    int periods[3] = {0, 0, 0};
    CPL::Cart_create (lammps->world, 3, npxyz_md, periods, myCoords.data(),
                      &icomm_grid);
   
    //TODO: get the origin from LAMMPS 

    double *xyz_orig = lammps->domain->boxlo;
    double *global_domain = lammps->domain->boxhi;
    //double xyz_orig[3] = {0.0 ,0.0, 0.0};

    //NOTE: Make sure set_timing is called before setup_cfd due to a unfixed bug
    //CPL::set_timing(initialstep, 0, 1.0);

    //Setup CPL run
    CPL::setup_md (icomm_grid, globaldomain, xyz_orig);

    // Store values of cell topology in socket
    getCellTopology();

    //Get timestep ratio and number of steps
    timestep_ratio = CPL::get<int> ("timestep_ratio");

}

void CPLSocketLAMMPS::getCellTopology() {

    // Cell sizes
    dx = CPL::get<double> ("dx");
    dy = CPL::get<double> ("dy");
    dz = CPL::get<double> ("dz");
       
    // Cell bounds for velocity BCs region
    CPL::get_bnry_limits(velBCRegion.data());
    CPL::my_proc_portion (velBCRegion.data(), velBCPortion.data());
    CPL::get_no_cells(velBCPortion.data(), velBCCells);

    // Cell bounds for the constrained region
    CPL::get_cnst_limits(cnstFRegion.data());
    CPL::my_proc_portion (cnstFRegion.data(), cnstFPortion.data());
    CPL::get_no_cells(cnstFPortion.data(), cnstFCells);

}




//Pack general using bitflag
void CPLSocketLAMMPS::allocateBuffers(const LAMMPS_NS::LAMMPS *lammps, int sendbitflag) {

    //Check what is to be packed and sent
    int packsize=0;
    if ((sendbitflag & VEL) == VEL){
        packsize += VELSIZE;
    }
    if ((sendbitflag & NBIN) == NBIN){
        packsize += NBINSIZE;
    }
    if ((sendbitflag & STRESS) == STRESS){
        packsize += STRESSSIZE;
        lammps->error->all(FLERR," sendbitflag stress not developed. Aborting.");
    }
    if ((sendbitflag & FORCE) == FORCE){
        packsize += FORCESIZE;
    }
    if ((sendbitflag & FORCECOEFF) == FORCECOEFF){
        packsize += FORCECOEFFSIZE;
    }
    if ((sendbitflag & VOIDRATIO) == VOIDRATIO){
        packsize += VOIDRATIOSIZE;
    }

    // LAMMPS computed velocity field
    int sendShape[4] = {packsize, velBCCells[0], velBCCells[1], velBCCells[2]};
    sendBuf.resize(4, sendShape);

    if (sendbitflag > 63)
        lammps->error->all(FLERR," sendbitflag bit flag unknown type. Aborting.");
}

void CPLSocketLAMMPS::setBndryAvgMode(int mode) {
	if (mode == AVG_MODE_ABOVE) {
		std::cout << "MODE ABOVE" << std::endl;
		bndry_shift_above = dy;
		bndry_shift_below = 0.0;
	}
	else if (mode == AVG_MODE_BELOW) {
		std::cout << "MODE BELOW" << std::endl;
		bndry_shift_above = 0.0;
		bndry_shift_below = dy;
	}
	else if (mode == AVG_MODE_MIDPLANE) {
		std::cout << "MODE MIDPLANE" << std::endl;
		bndry_shift_above = dy/2.0;
		bndry_shift_below = dy/2.0;
	}
	else if (mode == AVG_MODE_NONE) {
		std::cout << "MODE NO SHIFT" << std::endl;
		bndry_shift_above = 0.0;
		bndry_shift_below = 0.0;
	} else {
        lmp->error->all(FLERR,"Error - BndryAvg mode not specified");
    }
}

void CPLSocketLAMMPS::setupFixMDtoCFD(LAMMPS_NS::LAMMPS *lammps, int sendbitflag,
                                      int Nfreq, int Nrepeat, int Nevery)
{

    //Allocate buffers
    allocateBuffers(lammps, sendbitflag);

    double botLeft[3];
    CPL::map_cell2coord(velBCRegion[0] , velBCRegion[2], velBCRegion[4], botLeft);
    botLeft[1] -= bndry_shift_below;

    double topRight[3];
    CPL::map_cell2coord(velBCRegion[1] , velBCRegion[3], velBCRegion[5], topRight);


//    std::cout << "Extents " << botLeft[0] << " " <<  botLeft[1] << " " <<  botLeft[2] 
//                    << " " << topRight[0] << " " << topRight[1] << " " << topRight[2] << std::endl;

    topRight[0] += dx;
    topRight[1] += bndry_shift_above;
    topRight[2] += dz;

    double *global_domain = lammps->domain->prd;
    //double *global_domain = lammps->domain->boxhi;
//    std::cout << "Domain x = " << global_domain[0] << " CPL region x = " << topRight[0]-botLeft[0]
//              << " Domain y = " << global_domain[1] << " CPL region y = " << topRight[1]-botLeft[1] 
//              << " Domain z = " << global_domain[2] << " CPL region z = " << topRight[2]-botLeft[2] << std::endl;



//    std::cout << "Domain lo = " << lammps->domain->boxlo[0] << " "  
//                                << lammps->domain->boxlo[1] << " " 
//                                << lammps->domain->boxlo[2] << " "
//              << "Domain hi = " << lammps->domain->boxhi[0] << " "  
//                                << lammps->domain->boxhi[1] << " " 
//                                << lammps->domain->boxhi[2] << std::endl;

    //////////////////////////////////////////
    //This is the code sets the region
    //////////////////////////////////////////
    int ret;
    char topRight0str[20], topRight1str[20], topRight2str[20];
    char botLeft0str[20], botLeft1str[20], botLeft2str[20];
    ret = sprintf(topRight0str, "%f", topRight[0]);
    ret = sprintf(topRight1str, "%f", topRight[1]);
    ret = sprintf(topRight2str, "%f", topRight[2]);
    ret = sprintf(botLeft0str, "%f", botLeft[0]);
    ret = sprintf(botLeft1str, "%f", botLeft[1]);
    ret = sprintf(botLeft2str, "%f", botLeft[2]);

    // CFD BC region 
    char **regionarg = new char*[10];
    regionarg[0] = (char *) "cfdbcregion";
    regionarg[1] = (char *) "block";
    regionarg[2] = (char *) botLeft0str;
    regionarg[3] = (char *) topRight0str;
    regionarg[4] = (char *) botLeft1str;
    regionarg[5] = (char *) topRight1str;
    regionarg[6] = (char *) botLeft2str;
    regionarg[7] = (char *) topRight2str;
    regionarg[8] = (char *) "units";
    regionarg[9] = (char *) "box";
    lammps->domain->add_region(10, regionarg);
    delete [] regionarg;

    // CFD BC region 
    int iregion = lammps->domain->find_region("cfdbcregion");
    if (iregion < 0) lammps->error->all(FLERR,"Fix ID for iregion cfdbcregion does not exist");
    cfdbcregion = lammps->domain->regions[iregion];

    //////////////////////////////////////////
    //This code sets the compute
    //////////////////////////////////////////


    //New way using LAMMPS One
    std::stringstream cmp_str;
    cmp_str << std::setprecision(16) << "compute "  << "cfdbccompute "\
            << "all " << "chunk/atom bin/3d "\
            << "x lower " << dx << " "\
            << "y lower " << dy << " "\
            << "z lower " << dz << " "\
            << "ids every region " << "cfdbcregion "\
            << "units box " << " "\
            << "bound x " << botLeft[0] << " " << topRight[0] << " "\
            << "bound y " << botLeft[1] << " " << topRight[1] << " "\
            << "bound z " << botLeft[2] << " " << topRight[2];
    std::cout << "compute: " << cmp_str.str() << std::endl;
    lammps->input->one(cmp_str.str().c_str());

    //Old way passing a char array
    char dxstr[20], dystr[20], dzstr[20];
    char low_x[20], hi_x[20], low_y[20], hi_y[20], low_z[20], hi_z[20];
    ret = sprintf(dxstr, "%f", dx);
    ret = sprintf(dystr, "%f", dy);
    ret = sprintf(dzstr, "%f", dz);
    ret = sprintf(low_x, "%f", botLeft[0]);
    ret = sprintf(hi_x, "%f", topRight[0]);
    ret = sprintf(low_y, "%f", botLeft[1]);
    ret = sprintf(hi_y, "%f", topRight[1]);
    ret = sprintf(low_z, "%f", botLeft[2]);
    ret = sprintf(hi_z, "%f", topRight[2]);

    // CFD BC compute chunk 3d bins in y slice
    char **computearg = new char*[31];
    computearg[0] = (char *) "cfdbccompute";
    computearg[1] = (char *) "all";
    computearg[2] = (char *) "chunk/atom";
    computearg[3] = (char *) "bin/3d";
    computearg[4] = (char *) "x";
    computearg[5] = (char *) "lower";
    computearg[6] = (char *) dxstr;
    computearg[7] = (char *) "y";
    computearg[8] = (char *) "lower";
    computearg[9] = (char *) dystr;
    computearg[10] = (char *) "z";
    computearg[11] = (char *) "lower";
    computearg[12] = (char *) dzstr;
    computearg[13] = (char *) "ids";
    computearg[14] = (char *) "every";
    computearg[15] = (char *) "region";
    computearg[16] = (char *) "cfdbcregion";
    computearg[17] = (char *) "units";
    computearg[18] = (char *) "box";
    computearg[19] = (char *) "bound";
    computearg[20] = (char *) "x";
    computearg[21] = (char *) low_x;
    computearg[22] = (char *) hi_x;
    computearg[23] = (char *) "bound";
    computearg[24] = (char *) "y";
    computearg[25] = (char *) low_y;
    computearg[26] = (char *) hi_y;
    computearg[27] = (char *) "bound";
    computearg[28] = (char *) "z";
    computearg[29] = (char *) low_z;
    computearg[30] = (char *) hi_z;
    //lammps->modify->add_compute(31, computearg);
    //Get handle for compute
    int icompute = lammps->modify->find_compute("cfdbccompute");
    if (icompute < 0)
		lammps->error->all(FLERR,"Compute ID for compute cfdbccompute does not exist");
    cfdbccompute = lammps->modify->compute[icompute];
    std::cout << "cfdbccompute: " << cfdbccompute << " " << std::endl;
    
    delete [] computearg;

    //////////////////////////////////////////
    //This code sets the fix
    //////////////////////////////////////////
    // CFD BC averaging fix 
    // Average values are generated every Nfreq time steps, taken
    // from the average of the Nrepeat preceeding timesteps separated
    // by Nevery. For example, consider 
    // 
    //       Nfreq = 100;   Nrepeat = 6;   Nevery = 2;
    //
    // The average here would be taken over instantaneous snapshots at
    // timesteps 90, 92, 94, 96, 98, and 100 for the first value. The
    // next output would be an average from 190, 192, 194, 196, 198 and
    // 200 (and so on).

    //int Nfreq = timestep_ratio;
    //int Nrepeat = timestep_ratio;
    //int Nevery = 1;

    char Neverystr[20], Nrepeatstr[20], Nfreqstr[20];
    ret = sprintf(Neverystr, "%d", Nevery);
    ret = sprintf(Nrepeatstr, "%d", Nrepeat);
    ret = sprintf(Nfreqstr, "%d", Nfreq);

    std::stringstream fix_str;
    fix_str << "fix "  << "cfdbcfix "\
         << "all " << "ave/chunk "\
         << Neverystr << " " << Nrepeatstr << " "\
         << Nfreqstr << " "\
         << "cfdbccompute vx vy vz norm all";
    std::cout << "CPL: " << fix_str.str() << std::endl;
    lammps->input->one(fix_str.str().c_str());

    char **fixarg = new char*[14];
    fixarg[0] = (char *) "cfdbcfix";
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "ave/chunk";
    fixarg[3] = (char *) Neverystr;
    fixarg[4] = (char *) Nrepeatstr; 
    fixarg[5] = (char *) Nfreqstr; 
    fixarg[6] = (char *) "cfdbccompute";
    fixarg[7] = (char *) "vx";
    fixarg[8] = (char *) "vy";
    fixarg[9] = (char *) "vz";
    fixarg[10] = (char *) "norm";
    fixarg[11] = (char *) "all";
    fixarg[12] = (char *) "file";
    fixarg[13] = (char *) "cplchunk";
    lammps->modify->add_fix(14, fixarg);
    delete [] fixarg;

    //~ Set pointers for this newly-created fix
    int ifix = lammps->modify->find_fix("cfdbcfix");
    if (ifix < 0) lammps->error->all(FLERR,"Fix ID for fix cfdbcfix does not exist");
    cfdbcfix = lammps->modify->fix[ifix];
    
    //Move CPLSteps here to prevent a std::logic_error what():  basic_string::_S_construct null not valid
	int nsteps_md = CPL::get<int> ("nsteps_coupled") * timestep_ratio;
	std::string cmd =  "variable CPLSTEPS equal " + std::to_string(nsteps_md);
    lammps->input->one(cmd.c_str());
};

void CPLSocketLAMMPS::setupFixCFDtoMD(LAMMPS_NS::LAMMPS *lammps, 
                                      std::shared_ptr<std::string> forcetype, 
                                      std::vector<std::shared_ptr<std::string>> forcetype_args) {

    double botLeft[3];
    CPL::map_cell2coord(cnstFRegion[0] , cnstFRegion[2], cnstFRegion[4], botLeft);

    double topRight[3];
    CPL::map_cell2coord(cnstFRegion[1], cnstFRegion[3], cnstFRegion[5], topRight);
    topRight[0] += dx;
    topRight[1] += dy;
    topRight[2] += dz;

    // Tell LAMMPS to keep track of atoms in constrained region
    int ret;
    char topRight0str[20], topRight1str[20], topRight2str[20];
    char botLeft0str[20], botLeft1str[20], botLeft2str[20];
    ret = sprintf(topRight0str, "%f", topRight[0]);
    ret = sprintf(topRight1str, "%f", topRight[1]);
    ret = sprintf(topRight2str, "%f", topRight[2]);
    ret = sprintf(botLeft0str, "%f", botLeft[0]);
    ret = sprintf(botLeft1str, "%f", botLeft[1]);
    ret = sprintf(botLeft2str, "%f", botLeft[2]);

    char **regionarg = new char*[10];
    regionarg[0] = (char *) "cplforceregion";
    regionarg[1] = (char *) "block";
    regionarg[2] = (char *) botLeft0str;
    regionarg[3] = (char *) topRight0str;
    regionarg[4] = (char *) botLeft1str;
    regionarg[5] = (char *) topRight1str;
    regionarg[6] = (char *) botLeft2str;
    regionarg[7] = (char *) topRight2str;
    regionarg[8] = (char *) "units";
    regionarg[9] = (char *) "box";
    lammps->domain->add_region(10, regionarg);
    delete [] regionarg;

    int iregion = lammps->domain->find_region("cplforceregion");
    if (iregion < 0) lammps->error->all(FLERR,"Fix ID for iregion cplforceregion does not exist");
    cplforceregion = lammps->domain->regions[iregion];

    // CFD BC compute chunk 3d bins in y slice
    std::string cmd = "group cplforcegroup dynamic all region cplforceregion every 1";
    lammps->input->one (cmd.c_str());

    // Create a FixCPLForce instance
    std::string str = *forcetype;
    char * writable = new char[str.size() + 1];
    std::copy(str.begin(), str.end(), writable);
    writable[str.size()] = '\0'; // terminating 0
    char **fixarg = new char*[7+forcetype_args.size()];
    fixarg[0] = (char *) "cplforcefix";
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "cpl/force";
    fixarg[3] = (char *) "region";
    fixarg[4] = (char *) "cplforceregion"; 
    fixarg[5] = (char *) "forcetype";
    fixarg[6] = writable;
    int i=0;
    for(int i=0; i != forcetype_args.size(); i++) {
        auto arg = forcetype_args[i];
        std::string str(*arg);
        char * writable = new char[str.size() + 1];
        std::copy(str.begin(), str.end(), writable);
        writable[str.size()] = '\0'; // terminating 0
        fixarg[7+i] = writable;
    }
    lammps->modify->add_fix(7+forcetype_args.size(), fixarg);
    delete writable;
    delete [] fixarg;

    int ifix = lammps->modify->find_fix("cplforcefix");
    if (ifix < 0) lammps->error->all(FLERR,"Fix ID for fix cplforcefix does not exist");

    //Upcast Fix to child class FixCPLForce
    cplfix = dynamic_cast<FixCPLForce*>(lammps->modify->fix[ifix]);

    // Allocate received Buf field
    std::string fxyzType(*forcetype);
    int nval;
    if (fxyzType.compare("Flekkoy") == 0) {nval = 9; } 
    else if (fxyzType.compare("test") == 0) {nval = 3; } 
    else if (fxyzType.compare("Velocity") == 0) {nval = 3; } 
    else if (fxyzType.compare("Drag") == 0) {nval = 9; } 
    else if (fxyzType.compare("Di_Felice") == 0){nval = 9; } 
    else if (fxyzType.compare("Ergun") == 0) {nval = 9; }
    else if (fxyzType.compare("BVK") == 0) {nval = 9; }
    else {
        std::string cmd("CPLForce type ");
        cmd += fxyzType + " not defined";
        throw std::runtime_error(cmd);
    }
    int recvShape[4] = {nval, cnstFCells[0], cnstFCells[1], cnstFCells[2]};
    recvBuf.resize(4, recvShape);

	//Setup pointer to recieve buffer
	cplfix->setupBuf(recvBuf, cnstFPortion);

    //Again, no idea why this isn't automatically called by LAMMPS
    int vflag = 0;
    cplfix->setup(vflag);

}


//Pack general using bitflag
void CPLSocketLAMMPS::pack(const LAMMPS_NS::LAMMPS *lammps, int sendbitflag) {

        int *npxyz_md = lammps->comm->procgrid;
	    int nc_velBCRegion[3];
        CPL::get_no_cells(velBCRegion.data(), nc_velBCRegion);
        int row;
	    int glob_cell[3], loc_cell[3];
        double Vcell = dx*dy*dz;

        //Chosen arbitarily for now
        for (int i = velBCPortion[0]; i <= velBCPortion[1]; i++) {
        for (int j = velBCPortion[2]; j <= velBCPortion[3]; j++) {
	    for (int k = velBCPortion[4]; k <= velBCPortion[5]; k++) {  

		    glob_cell[0] = i; glob_cell[1] = j; glob_cell[2] = k;
		    CPL::map_glob2loc_cell(velBCPortion.data(), glob_cell, loc_cell);
            int ic=loc_cell[0]; int jc=loc_cell[1]; int kc=loc_cell[2];
            row = i*nc_velBCRegion[1]*nc_velBCRegion[2] + j*nc_velBCRegion[2] + k;
            int npack = 0;

            //Check what is to be packed and sent
            if ((sendbitflag & VEL) == VEL){
                //Get FSums internal to CPLForceTest
                std::string name("vSums");
                auto field_ptr = cplfix->fxyz->get_internal_fields(name);

                std::cout << "pack " <<  (field_ptr != nullptr) << " " << cfdbcfix->compute_array(row, 4) << " "  <<  std::endl;
                if (field_ptr != nullptr){
                    sendBuf(npack+0, ic, jc, kc) = field_ptr->get_array_value(0, ic, jc, kc);
                    sendBuf(npack+1, ic, jc, kc) = field_ptr->get_array_value(1, ic, jc, kc);
                    sendBuf(npack+2, ic, jc, kc) = field_ptr->get_array_value(2, ic, jc, kc);
                 } else {
                    double vx = cfdbcfix->compute_array(row, 4);  
                    double vy = cfdbcfix->compute_array(row, 5);  
                    double vz = cfdbcfix->compute_array(row, 6);

                    sendBuf(npack+0, ic, jc, kc) = vx;
                    sendBuf(npack+1, ic, jc, kc) = vy;
                    sendBuf(npack+2, ic, jc, kc) = vz; 
                    //lammps->error->all(FLERR," Array value vSums required by sendtype not collected in forcetype");
                }
                npack += VELSIZE;

                if (ic ==3 && jc == 3 && kc == 3){
                    std::cout << "CPLSocketLAMMPS::pack " << ic << " " << jc << " " << kc 
                              << " " <<  sendBuf(npack-3, ic, jc, kc)
                              << " " <<  sendBuf(npack-2, ic, jc, kc)
                              << " " <<  sendBuf(npack-1, ic, jc, kc) 
                              << " " << (field_ptr != nullptr) <<  std::endl;
                }

            }
            if ((sendbitflag & NBIN) == NBIN){
                //Get FSums internal to CPLForceTest
                std::string name("nSums");
                auto field_ptr = cplfix->fxyz->get_internal_fields(name);
                if (field_ptr != nullptr){
                    sendBuf(npack+0, ic, jc, kc) = field_ptr->get_array_value(0, ic, jc, kc);
                 } else {
                    double ncount = cfdbcfix->compute_array(row, 3);
                    sendBuf(npack, ic, jc, kc) = ncount; 
//                    lammps->error->all(FLERR," Array value nSums required by sendtype not collected in forcetype");
                }
                npack += NBINSIZE;
            }
            if ((sendbitflag & STRESS) == STRESS){
                lammps->error->all(FLERR," sendbitflag stress not developed. Aborting.");
                npack += STRESSSIZE;
            }
            if ((sendbitflag & FORCE) == FORCE){

                //Get FSums internal to CPLForceTest
                std::string name("FSums");
                auto field_ptr = cplfix->fxyz->get_internal_fields(name);
                if (field_ptr != nullptr){
                    sendBuf(npack+0, ic, jc, kc) = field_ptr->get_array_value(0, ic, jc, kc);
                    sendBuf(npack+1, ic, jc, kc) = field_ptr->get_array_value(1, ic, jc, kc);
                    sendBuf(npack+2, ic, jc, kc) = field_ptr->get_array_value(2, ic, jc, kc);
                 } else {
                    lammps->error->all(FLERR," Array value FSums required by sendtype not collected in forcetype");
                }
                npack += FORCESIZE;
            }
            if ((sendbitflag & FORCECOEFF) == FORCECOEFF){

                std::string name("FcoeffSums");
                auto field_ptr = cplfix->fxyz->get_internal_fields(name);
                if (field_ptr != nullptr){
                    sendBuf(npack+0, ic, jc, kc) = field_ptr->get_array_value(0, ic, jc, kc);
                } else {
                    lammps->error->all(FLERR," Array value FcoeffSums required by sendtype not collected in forcetype");
                }
                npack += FORCECOEFFSIZE;
            }
            if ((sendbitflag & VOIDRATIO) == VOIDRATIO){

                std::string name("volSums");
                auto field_ptr = cplfix->fxyz->get_internal_fields(name);
                if (field_ptr != nullptr){
                    //Send sum of volume directly
                    sendBuf(npack, ic, jc, kc) = field_ptr->get_array_value(0, ic, jc, kc);
//                    double phi = field_ptr->get_array_value(0, ic, jc, kc)/Vcell;
//                    if (phi > 1.) {
//                        //std::cout << "Warning, eps = 0 so set to 0.1 in CPLSocketLAMMPS::packGran" << std::endl;
//                        sendBuf(npack, ic, jc, kc) = 0.01;
//                    } else {
//                        sendBuf(npack, ic, jc, kc) = 1.0 - phi;
//                        //if (phi != 0.0)
//                        //    std::cout << "CPLSocketLAMMPS::pack " << ic << " " << jc << " " << kc << " " <<  sendBuf(npack, ic, jc, kc) << std::endl;
//                    }
                } else {
                    lammps->error->all(FLERR," Array value volSums required by sendtype not collected in forcetype");
                }
                //std::cout << i << " " << j << " " << k << " " << 1. - Granfxyz.volSums(i,j,k)/Vcell << std::endl;
                npack += VOIDRATIOSIZE;
            }
        }}}

}
    
    
void CPLSocketLAMMPS::send() {

    // Send the data to CFD
//    std::cout << "CPLSocketLAMMPS::send "
//              << " " << sendBuf.shape(0) 
//              << " " << sendBuf.shape(1) 
//              << " " << sendBuf.shape(2) 
//              << " " << sendBuf.shape(3) << std::endl;
    CPL::send(sendBuf.data(), sendBuf.shapeData(), velBCRegion.data());
};

void CPLSocketLAMMPS::receive() {
    // Receive from CFD
    CPL::recv(recvBuf.data(), recvBuf.shapeData(), cnstFRegion.data());
//    std::cout << "CPLSocketLAMMPS::recv "
//              << " " << recvBuf.shape(0) 
//              << " " << recvBuf.shape(1) 
//              << " " << recvBuf.shape(2) 
//              << " " << recvBuf.shape(3) << std::endl;

};


////General function to parse args and return a vector...

//std::vector<std::shared_ptr<std::string>> parse_arguments(int narg, char **arg, 
//                                                            std::vector<std::string> endargs) {

//    std::shared_ptr<std::string> forcetype;
//    std::vector<std::shared_ptr<std::string>> forcetype_args;

//    for (int iarg=0; iarg<narg; iarg+=1){
//        std::string arguments(arg[iarg]);
//        if (arguments == "forcetype"){
//            if (iarg+1<narg) {
//                for (int jarg=iarg+1; jarg<narg; jarg+=1){
//                    std::shared_ptr<std::string> forcetype_arg;
//                    forcetype_arg = std::make_shared<std::string>(arg[jarg]);
//                    //Check if we have read another argument type
//                    std::string forceType_arg(*forcetype_arg);
//                    for ( auto &endarg : endargs ) {
//                        if (  forceType_arg.compare(endarg) == 0)
//                            break;
//                    }
//                    //Otherwise it is a sendtype argument and should be added
//                    std::string forceType(*forcetype);
//                    std::cout << "Lammps FixCPLForce forcetype: "  << forceType << " with args "
//                              <<  forceType_arg << std::endl;
//                    forcetype_args.push_back(forcetype_arg);
//                }
//                
//            }
//        }
//    }
//    return forcetype_args;
//}






