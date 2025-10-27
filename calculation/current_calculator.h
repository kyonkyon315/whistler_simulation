#ifndef CURRENT_CALCULATOR_H
#define CURRENT_CALCULATOR_H
#include "distribution_function.h"
#include "current.h"
#include "parameters.h"
class CurrentCalculator{
    using Value = double;
    using Index = int;

private:
    const AxisInfo_x& axis_x;
    const AxisInfo_r& axis_r;
    const AxisInfo_theta& axis_t;
    const AxisInfo_phi& axis_p;

public:
    CurrentCalculator(
        const AxisInfo_x& axis_x,
        const AxisInfo_r& axis_r,
        const AxisInfo_theta& axis_t,
        const AxisInfo_phi& axis_p
    ):
        axis_x(axis_x),
        axis_r(axis_r),
        axis_t(axis_t),
        axis_p(axis_p)
    {}
    //fをもとに電流を計算し、targetに格納する。
    void calc(const DistributionFunction& f,Current& target){
        Index x_size = axis_x.get_num_grid();
        Index r_size = axis_r.get_num_grid();
        Index t_size = axis_t.get_num_grid();
        Index p_size = axis_p.get_num_grid();
        
        #pragma omp parallel for
        for(Index x=0;x<x_size;++x){
            target.at<Direction::x>(x) = 0.;
            target.at<Direction::y>(x) = 0.;
            target.at<Direction::z>(x) = 0.;

            for(Index r=0;r<r_size;++r){
                for(Index t=0;t<t_size;++t){
                    for(Index p=0;p<p_size;++p){
                        target.at<Direction::x>(x)
                            +=  f.dV_at(x,r,t,p)*
                                f.get_velocity_at<Direction::x>(x,r,t,p);
                        target.at<Direction::y>(x)
                            +=  f.dV_at(x,r,t,p)*
                                f.get_velocity_at<Direction::y>(x,r,t,p);
                        target.at<Direction::z>(x)
                            +=  f.dV_at(x,r,t,p)*
                                f.get_velocity_at<Direction::z>(x,r,t,p);
                    }
                }
            }
            Value q = Parameters::elementary_charge;
            target.at<Direction::x>(x)*=q;
            target.at<Direction::y>(x)*=q;
            target.at<Direction::z>(x)*=q;
        }
    }
};
#endif// CURRENT_CALCULATOR_H
