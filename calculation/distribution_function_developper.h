#ifndef DISTRIBUTION_FUNCTION_DEVELOPPER_H
#define DISTRIBUTION_FUNCTION_DEVELOPPER_H
#include "distribution_function.h"
#include "current.h"
#include "electro_magnetic_field.h"
#include <string>
class DistributionFunctionDevelopper{
    public:

    //df/dtをtargetに格納する。
    void deffirenciate(
        const DistributionFunction& f,
        const ElectroMagneticField& field,
        const Current& current,
        const std::string& type/*"main" or "half"*/,
        DistributionFunction& target
    ){

    }

};

#endif// DISTRIBUTION_FUNCTION_DEVELOPPER_H
