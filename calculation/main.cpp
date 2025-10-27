#include <iostream>
#include "axis_info.h"
#include "distribution_function.h"
#include "electro_magnetic_field.h"
#include "current.h"
#include "developper.h"
using Value = double;
using Index = int;

int main(){
    Index x_num_grid=1000;
    Index vr_num_grid=100;
    Index vtheta_num_grid=100;
    Index vphi_num_grid=100;
    Value x_begin=0.,x_end=100.;
    Value vr_end=100.;

    AxisInfo_x x(x_begin,x_end,x_num_grid);
    AxisInfo_r v_r(vr_end,vr_num_grid);
    AxisInfo_theta v_theta(vtheta_num_grid);
    AxisInfo_phi v_phi(vphi_num_grid);


    DistributionFunction f(x,v_r,v_theta,v_phi);
    ElectroMagneticField field(x);
    std::cout<<v_phi[50]<<std::endl;

    Current current(x);
    std::cout<<v_phi[50]<<std::endl;

    current.at<Direction::x>(4);

    std::cout<<field.at<FieldType::B,Direction::x>(4)<<std::endl;


    Developper developper(f,x,v_r,v_theta,v_phi);

    Index num_time_step=1000;
    Value dt=0.1;
    for(Index step=0;step<num_time_step;++step){
        std::cout<<step<<"\n";
        developper.update(f,field,current,dt);
    }

    return 0;
}
