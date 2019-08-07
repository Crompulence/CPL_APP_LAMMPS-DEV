#ifdef FIX_CLASS

FixStyle(cpl/bc, FixCPLBc)

#else

#ifndef LMP_FIX_CPL_BC_H
#define LMP_FIX_CPL_BC_H

#include "fix.h"
#include "cpl/cpl.h"
#include "CPLSocketLAMMPS.h"
#include <memory>
#include "LAMMPSDep.h"
#include "fix_cpl_init.h"

class FixCPLBc: public LAMMPS_NS::Fix {

public:

    FixCPLBc(class LAMMPS_NS::LAMMPS *lammps, int narg, char **arg);
	~FixCPLBc();
    void post_constructor();
    int setmask(){return 0;}

    // Pool of boundary conditions for CFD
    CPL::OutgoingFieldPool* bcPool;
    // This points to FixCPLInit depPool object.
    DepPoolT* depPool;

private:
    FixCPLInit* fixCPLInit;

};

DEPFUNC_DEF(cfdbcregion_depfunc);
DEPFUNC_DEF(cfdbc_chunks_depfunc);
DEPFUNC_DEF(cfdbc_property_depfunc);
DEPFUNC_DEF(cfdbc_vcom_depfunc);
DEPFUNC_DEF(cfdbc_tpartial_depfunc);
DEPFUNC_DEF(cfdbc_temp_depfunc);
DEPFUNC_DEF(cfdbc_fix_depfunc);

#endif
#endif
