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
    
    DistributionFunctionDevelopper distribution_function_developper;
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
        distribution_function_developper(
            axis_info_x,
            axis_info_vr,
            axis_info_v_theta,
            axis_info_v_phi),
        current_calculator(
            axis_info_x,
            axis_info_vr,
            axis_info_v_theta,
            axis_info_v_phi),
        fdtd_solver()
    {
    }

   

    void update(
        DistributionFunction& f,
        ElectroMagneticField& field,
        Current& current,
        Value dt
    ){
    /*
    splitting 法により計算
    # 1: x 半ステップ
    f = Advect_x(f, dt/2)

    # 2: compute J (optionally), advance fields half step
    J = compute_current(f)
    E,B = FDTD_update(E,B,J, dt/2)

    # 3: E 半ステップ
    f = Advect_v_by_E(f, E, dt/2)

    # 4: B フルステップ
    f = Advect_v_by_B(f, B, dt)   # 速度空間での回転（逆追跡＋補間）

    # 5: E 半ステップ
    f = Advect_v_by_E(f, E, dt/2)

    # 6: compute J, advance fields last half step
    J = compute_current(f)
    E,B = FDTD_update(E,B,J, dt/2)

    # 7: x 半ステップ
    f = Advect_x(f, dt/2)

    # now f,E,B are at t^{n+1}
    
    */
        distribution_function_developper.advect_x(f,dt/2);
        //fをv・∇fの寄与に基づきdt/2発展

        current_calculator.calc(f,current);
        //電流を計算してcurrentに格納

        fdtd_solver.develop(f,current,dt/2,field);
        //fieldをfdtdに基づきdt/2 発展

        distribution_function_developper.advect_E(f,field,dt/2);
        //fをqE・∇_v fの寄与に基づきdt/2発展

        distribution_function_developper.advect_B(f,field,dt);
        //fをqVxB・∇_v fの寄与に基づきdt発展
        
        distribution_function_developper.advect_E(f,field,dt/2);
        //fをqE・∇_v fの寄与に基づきdt/2発展

        current_calculator.calc(f,current);
        //電流を計算してcurrentに格納

        fdtd_solver.develop(f,current,dt/2,field);
        //fieldをfdtdに基づきdt/2 発展

        distribution_function_developper.advect_x(f,dt/2);
        //fをv・∇fの寄与に基づきdt/2発展
    }
};

#endif// DEVELOPPER_H
