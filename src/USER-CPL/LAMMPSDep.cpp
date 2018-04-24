#include "LAMMPSDep.h"
#include "modify.h"
#include "input.h"
#include "error.h"
#include "group.h"


void DepLoader::loadDeps(const DepListT& dep_list, DepPoolT& dep_pool) {
    DepListT::const_iterator i;
    for (i = dep_list.begin(); i != dep_list.end(); i++)
        dep_pool[*i]->load();
}


LAMMPSDep::LAMMPSDep(std::string dep_name, const DepListT& dep_list, 
                     void* obj, 
                     LAMMPS_NS::LAMMPS* lammps, DepFuncT dep_func,
                     int n_every) : name(dep_name), 
                     depList(dep_list), cmd(dep_func), nevery(n_every), 
                     lmp(lammps), loaded(false) {

    elem_name = dep_name;
    //LAMMPS Specific
    cmd_str = cmd(name, obj, lammps, nevery);
}


bool LAMMPSDep::load() { 
    if (!loaded) { 
        loadDeps(depList, *pool);
        //Specific to LAMMPS
        lmp->input->one(cmd_str.c_str());
        load_();
        loaded = true;
        return true;
    }
    return false;

}

    
void LAMMPSDepFix::load_() { 
    int ifix = lmp->modify->find_fix(name.c_str());
    if (ifix < 0) {
        std::string err_msg = "Fix ID for fix " + name +\
                               " does not exist. An error has ocurred.";
        lmp->error->all(FLERR,err_msg.c_str());
    }
    fix = lmp->modify->fix[ifix];
}


void LAMMPSDepGroup::load_() {
    groupId = lmp->group->find(name.c_str());
    if (groupId < 0) {
        std::string err_msg = "Group ID for group " + name +\
                               " does not exist. An error has ocurred.";
        lmp->error->all(FLERR,err_msg.c_str());
    }
}

void LAMMPSDepRegion::load_() { 
    // In domain.h, find_region(char*) should be find_region(const char*) I think.
    int iregion = lmp->domain->find_region(const_cast<char*>(name.c_str()));
    if (iregion < 0) {
        std::string err_msg = "Region ID for region " + name +\
                               " does not exist. An error has ocurred.";
        lmp->error->all(FLERR,err_msg.c_str());
    }
    region = lmp->domain->regions[iregion];
}

void LAMMPSDepCompute::load_() { 
    int icompute = lmp->modify->find_compute(name.c_str());
    if (icompute < 0) {
        std::string err_msg = "Compute ID for compute " + name +\
                               " does not exist. An error has ocurred.";
        lmp->error->all(FLERR,err_msg.c_str());
    }
    compute = lmp->modify->compute[icompute];
}


