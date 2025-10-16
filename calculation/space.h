#ifndef SPACE_H
#define SPACE_H
#include "axis_info.h"
#include <vector>
class Space{
    using Value = double;
    using Index = int;

    public:
    const AxisInfo_x& axis_x;
    const AxisInfo_r& axis_r;
    const AxisInfo_theta& axis_t;
    const AxisInfo_phi& axis_p;
    const Index overlap;
    const Index r_theta_phi;
    const Index theta_phi;
    const Index phi;
    std::vector<Value> dV;
    std::vector<Value> velocity;
    Space(
        const AxisInfo_x&       axis_x,
        const AxisInfo_r&       axis_r,
        const AxisInfo_theta&   axis_t,
        const AxisInfo_phi&     axis_p,
        Index overlap=5
    ):
        axis_x(axis_x),
        axis_r(axis_r),
        axis_t(axis_t),
        axis_p(axis_p),
        overlap(overlap),
        r_theta_phi(
            (axis_r.get_num_grid()+2*overlap)*
            (axis_t.get_num_grid()+2*overlap)*
            (axis_p.get_num_grid()+2*overlap)
        ),
        theta_phi(
            (axis_t.get_num_grid()+2*overlap)*
            (axis_p.get_num_grid()+2*overlap)
        ),
        phi(
            (axis_p.get_num_grid()+2*overlap)
        )
    {
        pre_calc_dV();
    }

    Index id2vec_id(Index x,Index r,Index t,Index p)const{
        
    }

    private:
    void pre_calc_dV(){
        for(Index r=0;r<axis_r.get_num_grid();++r){
            for(Index t=0;t<axis_t.get_num_grid();++t){
                for(Index p=0;p<axis_p.get_num_grid();++p){
                    
                }
            }
        }
    }
};
#endif// SPACE_H