#ifndef LAMMPS_INCOMING_FIELD_H
#define LAMMPS_INCOMING_FIELD_H

#include "cpl/CPL_field.h"
#include "LAMMPSDep.h"
#include "cpl/TransmittingField.h"
#include <valarray>

class LAMMPSIncomingField : public CPL::IncomingField, public DepLoader {
    public:
        LAMMPSIncomingField() : IncomingField() {}
        LAMMPSIncomingField(std::string name, const PortionField& portion_field,
                            const PortionField& field, const DepListT& dep_list,
                            DepPoolT& dep_pool, LAMMPS_NS::LAMMPS* lammps) :
                            CPL::IncomingField(name, portion_field, field), 
                            lmp(lammps), depPool(&dep_pool),
                            depList(dep_list) {}

        LAMMPS_NS::LAMMPS* lmp;
        DepListT depList;
        DepPoolT* depPool;
        virtual ~LAMMPSIncomingField(){}

};


class StressIncomingField:  public LAMMPSIncomingField {
    public:
        using LAMMPSIncomingField::LAMMPSIncomingField;
        void setup();
};


#endif //LAMMPS_INCOMING_FIELD_H
