#ifndef DEVELOPPER_H
#define DEVELOPPER_H
#include "axis_info.h"
#include "distribution_function.h"
#include "electro_magnetic_field.h"
#include "current.h"

#include "distribution_function_developper.h"
#include "current_calculator.h"
class Developper{

    using Value = double;
    using Index = int;

    private:
    DistributionFunction k1,k2;
    DistributionFunctionDevelopper developper;
    CurrentCalculator current_calculator;
    FDTDSolver  fdtd_solver;


    public:
    Developper(const DistributionFunction& f):
        k1(f),k2(f)
    {}

    void update(
        DistributionFunction& f,
        ElectroMagneticField& field,
        Current& current,
        Value dt
    ){
        developper.deffirenciate(f,field,current,"main",k1);//k1=df/dt
        k1.mul(dt/2.);//k1=k1*dt/2  
        k1.add(f);//k1=f+k1  k1=f(t+dt/2)
        current_calculator.calc(k1,current);//k1->current(t+dt/2)
        fdtd_solver.calc(k1,current,field,"main");//field(t+dt)を計算
        developper.deffirenciate(k1,field,current,"half",k2);//k2=dfdt(t+dt/2)
        k2.mul(dt);//k2=k2*dt
        f.add(k2);
        fdtd_solver.calc(f,current,field,"half")//current(t+dt)を計算

    }
};

#endif// DEVELOPPER_H
