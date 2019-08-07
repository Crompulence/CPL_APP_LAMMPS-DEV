#ifndef LAMMPS_OUTGOING_FIELD_H
#define LAMMPS_OUTGOING_FIELD_H

#include "cpl/CPL_field.h"
#include "LAMMPSDep.h"
#include "cpl/TransmittingField.h"
#include <valarray>

class LAMMPSOutgoingField : public CPL::OutgoingField, public DepLoader {
    public:
        LAMMPSOutgoingField() : OutgoingField() {}
        LAMMPSOutgoingField(std::string name,
                            const DepListT& dep_list, DepPoolT* dep_pool, 
                            LAMMPS_NS::LAMMPS* lammps) :
                            CPL::OutgoingField(name), lmp(lammps),
                            depPool(dep_pool), depList(dep_list) {}

        LAMMPSOutgoingField(std::string name, const PortionField& portion_field,
                            const PortionField& field, const DepListT& dep_list,
                            DepPoolT* dep_pool, LAMMPS_NS::LAMMPS* lammps) :
                            CPL::OutgoingField(name, portion_field, field),
                            lmp(lammps), depPool(dep_pool), depList(dep_list) {}



        LAMMPS_NS::LAMMPS* lmp;
        DepListT depList;
        DepPoolT* depPool;
        std::string bc_dimension;
        virtual ~LAMMPSOutgoingField(){}

};


class VelOutgoingField:  public LAMMPSOutgoingField {
    public:
        using LAMMPSOutgoingField::LAMMPSOutgoingField;
        void pack_(const std::vector<int>& glob_cell, 
                   const std::vector<int>& loc_cell,
                   const std::valarray<double>& coord);
        void setup();
        void update();
        std::valarray<double> vAvg;
        bool averageVels;
};

class TemperatureOutgoingField:  public LAMMPSOutgoingField {
    public:
        using LAMMPSOutgoingField::LAMMPSOutgoingField;
        void pack_(const std::vector<int>& glob_cell, 
                   const std::vector<int>& loc_cell,
                   const std::valarray<double>& coord);
        void setup();
};


class NbinOutgoingField:  public LAMMPSOutgoingField {
    public:
        using LAMMPSOutgoingField::LAMMPSOutgoingField;
        void pack_(const std::vector<int>& glob_cell,
                   const std::vector<int>& loc_cell,
                   const std::valarray<double>& coord);
        void setup();
};




#endif //LAMMPS_OUTGOING_FIELD_H
