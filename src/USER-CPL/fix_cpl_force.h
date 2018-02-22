#ifdef FIX_CLASS

FixStyle(cpl/force, FixCPLForce)

#else

#ifndef LMP_FIX_CPL_FORCE_H
#define LMP_FIX_CPL_FORCE_H

#include<memory>

#include "fix.h"
//#include "fix_rigid.h"

#include "cpl/cpl.h"
#include "cpl/CPL_ndArray.h"
#include "cpl/CPL_force.h"

class FixCPLForce : public LAMMPS_NS::Fix {

public:

    FixCPLForce(class LAMMPS_NS::LAMMPS *lammps, int narg, char **arg);
    int setmask();
	void setup(int vflag);
	void setupBuf(CPL::ndArray<double>& Buf, std::vector<int>& portion);
    void apply();
    void updateProcPortion (std::vector<int>& portion);
    std::shared_ptr<std::string> forcetype;
    std::unique_ptr<CPLForce> fxyz;

private:

	CPL::ndArray<double>* cfdBuf;
    std::vector<int> procPortion;
    std::vector<double> fi;

    std::vector<std::shared_ptr<std::string>> forcetype_args;
    //std::unique_ptr<LAMMPS_NS::FixRigid> clumpfix;

};

#endif
#endif









// ****************************************************************
// Copied here for record of what not to do! Rigid mols are local!
// and rigid molecules are private so not sure how we get at them



    //If we have rigid bodies, need to include these too
    //int ifix = lmp->modify->find_fix("clumps");
    //if (ifix < 0) lmp->error->all(FLERR,"Fix ID for fix clumps does not exist");

    //auto clumpfix = dynamic_cast<LAMMPS_NS::FixRigid>(lmp->modify->fix[ifix]);
    //int nbody = lmp->modify->fix[ifix]->nbody;
    //int *nrigid = lmp->modify->fix[ifix]->nrigid; //Number of molecules in rigid body
    //int *body2mol = lmp->modify->fix[ifix]->body2mol;

        //LOOP over all rigid bodies and get forces
//        // This would not work as these are protected
//        for (int b = 0; b < nbody; ++b)
//        {
//            //Get number of molecules in rigid body
//        	for (int bi = 0; bi < nrigid[b]; ++bi)
//        	{

//                //Get global molecular id from body id
//                int i = body2mol[bi];
//           		if (mask[i] & groupbit)
//            	{

//		            //Get local molecule data
//		            mi = rmass[i];
//		            radi = radius[i];
//		            for (int n=0; n<3; n++){
//		                xi[n]=x[i][n]; 
//		                vi[n]=v[i][n]; 
//		                ai[n]=f[i][n];
//		            }

//		            // Sum all the weights for each cell.
//		            fxyz->pre_force(xi, vi, ai, mi, radi, pot);

//            	}
//            }
//        }
