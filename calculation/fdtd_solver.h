#ifndef FDTD_SOLVER_H
#define FDTD_SOLVER_H
#include "distribution_function.h"
#include "current.h"
#include "electro_magnetic_field.h"
class FDTDSolver{
public:
    void develop(
        const DistributionFunction& f,
        const Current& current,
        Value dt,
        ElectroMagneticField& electro_magnetic_field
    ){

    }
};
#endif// FDTD_SOLVER_H
