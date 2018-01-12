#ifdef FIX_CLASS

FixStyle(cpl/constrain, FixCPLConstrain)

#else

#ifndef LMP_FIX_CPL_CONSTRAIN_H
#define LMP_FIX_CPL_CONSTRAIN_H

#include "fix.h"
#include "cpl/cpl.h"
#include "CPLSocketLAMMPS.h"
#include <memory>
#include "LAMMPSDep.h"
#include "fix_cpl_init.h"

class FixCPLConstrain: public LAMMPS_NS::Fix {

public:

    FixCPLConstrain(class LAMMPS_NS::LAMMPS *lammps, int narg, char **arg);
	~FixCPLConstrain(){};
    void post_constructor();
    int setmask(){}

    // Pool of constrains for MD
    CPL::IncomingFieldPool* cnstPool;
    DepPoolT* depPool;

private:
    FixCPLInit* fixCPLInit;

};

DEPFUNC_DEF(cplforceregion_depfunc);
DEPFUNC_DEF(cplforcegroup_depfunc);
DEPFUNC_DEF(cplforcefix_depfunc);

#endif
#endif
