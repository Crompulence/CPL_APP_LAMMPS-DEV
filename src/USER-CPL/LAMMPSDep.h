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

    "LammpsDep" class which encapsulate a lammps group, region, fix, compute, etc.

Author(s)

   Eduardo Ramos 

*/
#ifndef CPL_LAMMPS_DEP_H_INCLUDED
#define CPL_LAMMPS_DEP_H_INCLUDED

#include <vector>
#include "mpi.h"
#include "lammps.h"
#include "cpl/cpl.h"
#include "region.h"
#include <string>
#include "CPLSocketLAMMPS.h"
#include "cpl/PoolElement.h"


#define DEPFUNC_DEF(fname) std::string fname(std::string dep_name, \
                          void* obj, \
                          LAMMPS_NS::LAMMPS* lmp, int nevery=-1)
#define DEPFUNC_IMP(fname) std::string fname(std::string dep_name, \
                           void* obj, \
                           LAMMPS_NS::LAMMPS* lmp, int nevery)

class LAMMPSDep;

typedef std::vector<std::string> DepListT;
typedef CPL::Pool<LAMMPSDep> DepPoolT;
typedef std::string (*DepFuncT)(std::string, void* obj, LAMMPS_NS::LAMMPS*, int);

class DepLoader {
    public:
        DepLoader(){}
        virtual ~DepLoader(){};
        void loadDeps(const DepListT& dep_list, DepPoolT& dep_pool);
};

class LAMMPSDep : public DepLoader, public CPL::PoolElement<LAMMPSDep>{

public:
    // Construct from no arguments
    LAMMPSDep(std::string dep_name, const DepListT& dep_list, 
              void* obj, 
              LAMMPS_NS::LAMMPS* lammps, DepFuncT dep_func, 
              int n_every = -1); 
    LAMMPSDep(){}
    std::string name;
    DepListT depList;
    int nevery;
    std::string depType;
    bool loaded;
    bool load();

    //specific to lammps
    DepFuncT cmd;
    std::string cmd_str;
    LAMMPS_NS::LAMMPS* lmp;
    virtual ~LAMMPSDep(){};

    // Virtual members
    virtual void load_()=0;
};

class LAMMPSDepFix : public LAMMPSDep {
    public:
        using LAMMPSDep::LAMMPSDep;
        LAMMPS_NS::Fix* fix;
        virtual void load_();
        virtual ~LAMMPSDepFix(){std::cout << "DEPFIX DESTRUCTOR" << std::endl;};
};

class LAMMPSDepGroup : public LAMMPSDep {
    public:
        using LAMMPSDep::LAMMPSDep;
        // In this case all groups are handled under a single class
        // and arrays like group->bitmask, group->dynamic, group->names
        // are accessed with the group id which we save.
        int groupId;
        virtual void load_();
        virtual ~LAMMPSDepGroup(){};
};

class LAMMPSDepRegion : public LAMMPSDep {
    public:
        using LAMMPSDep::LAMMPSDep;
        LAMMPS_NS::Region* region;
        virtual void load_();
        virtual ~LAMMPSDepRegion(){std::cout << "DEPREGION DESTRUCTOR" << std::endl;};
};

class LAMMPSDepCompute : public LAMMPSDep {
    public:
        using LAMMPSDep::LAMMPSDep;
        LAMMPS_NS::Compute* compute;
        virtual void load_();
        virtual ~LAMMPSDepCompute(){std::cout << "DEPCOMPUTE DESTRUCTOR" << std::endl;};
};





#endif //CPL_LAMMPS_DEP_H_INCLUDED
