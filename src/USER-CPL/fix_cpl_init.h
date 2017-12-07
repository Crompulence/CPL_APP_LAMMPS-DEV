#ifdef FIX_CLASS

FixStyle(cpl/init, fixCPLInit)

#else

#ifndef LMP_FIX_CPL_INIT_H
#define LMP_FIX_CPL_INIT_H

#include "fix.h"
#include "cpl/cpl.h"
#include "CPLSocketLAMMPS.h"
#include <memory>
#include "LAMMPSDep.h"

class fixCPLInit : public LAMMPS_NS::Fix {

public:

    fixCPLInit(class LAMMPS_NS::LAMMPS *lammps, int narg, char **arg);
	~fixCPLInit();
    int setmask();
    void setup (int vflag); 
    void end_of_step();
    void post_constructor();
    void setas_last_fix();
	void init();

    // Pool of boundary conditions for CFD
    CPL::OutgoingFieldPool bcPool;
    CPL::IncomingFieldPool cnstPool;
    
    std::string forceType;
    std::string sendType;
    std::string bndryAvg;
    DepPoolT depPool;

private:
    LAMMPS_NS::LAMMPS* lmp;
    void setupDeps();
    void setupBcs();


};

DEPFUNC_DEF(cfdbcregion_depfunc);
DEPFUNC_DEF(cfdbccompute_depfunc);
DEPFUNC_DEF(cfdbcfix_depfunc);
DEPFUNC_DEF(cplforceregion_depfunc);
DEPFUNC_DEF(cplforcegroup_depfunc);
DEPFUNC_DEF(cplforcefix_depfunc);

#endif
#endif
