#ifndef DEVELOPPER_H
#define DEVELOPPER_H
#include "axis_info.h"
#include "distribution_function.h"
#include "electro_magnetic_field.h"
#include "current.h"

#include "distribution_function_developper.h"
#include "current_calculator.h"
#include "fdtd_solver.h"
class Developper{

    using Value = double;
    using Index = int;

private:
    DistributionFunction k1;
    DistributionFunction k2;
    DistributionFunctionDevelopper developper;
    CurrentCalculator current_calculator;
    FDTDSolver  fdtd_solver;


public:
    Developper(
        const DistributionFunction& f,
        const AxisInfo_x& axis_info_x,
        const AxisInfo_r& axis_info_vr,
        const AxisInfo_theta& axis_info_v_theta,
        const AxisInfo_phi& axis_info_v_phi
    ): 
        k1(f),            // ← コピーコンストラクタを使って初期化
        k2(f),
        developper(),     // ← デフォルトコンストラクタを呼ぶ（必要なら引数付きに変更）
        current_calculator(axis_info_x,axis_info_vr,axis_info_v_theta,axis_info_v_phi),
        fdtd_solver()
    {
    }

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
