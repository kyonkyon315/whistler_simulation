#ifndef DISTRIBUTION_FUNCTION_H
#define DISTRIBUTION_FUNCTION_H
#include <vector>
#include "axis_info.h"
class DistributionFunction{

    using Value = double;
    using Index = int;

    private:
    const AxisInfo_x& axis_info_x;
    const AxisInfo_r& axis_info_vr;
    const AxisInfo_theta& axis_info_v_theta;
    const AxisInfo_phi& axis_info_v_phi;
    const Index r_theta_phi;
    const Index theta_phi;
    const Index phi;
    const Index overlap;
    std::vector<Value> f;
    std::vector<Value> dV_table;

    public:
    DistributionFunction(
        const AxisInfo_x&       axis_info_x,
        const AxisInfo_r&       axis_info_vr,
        const AxisInfo_theta&   axis_info_v_theta,
        const AxisInfo_phi&     axis_info_v_phi,
        Index overlap=5
    ):
        axis_info_x(axis_info_x),
        axis_info_vr(axis_info_vr),
        axis_info_v_theta(axis_info_v_theta),
        axis_info_v_phi(axis_info_v_phi),
        overlap(overlap),
        r_theta_phi(
            (axis_info_vr.get_num_grid()+2*overlap)*
            (axis_info_v_theta.get_num_grid()+2*overlap)*
            (axis_info_v_phi.get_num_grid()+2*overlap)
        ),
        theta_phi(
            (axis_info_v_theta.get_num_grid()+2*overlap)*
            (axis_info_v_phi.get_num_grid()+2*overlap)
        ),
        phi(
            (axis_info_v_phi.get_num_grid()+2*overlap)
        )
    {
        f.resize(
            (axis_info_x.get_num_grid()+2*overlap)*
            (axis_info_vr.get_num_grid()+2*overlap)*
            (axis_info_v_theta.get_num_grid()+2*overlap)*
            (axis_info_v_phi.get_num_grid()+2*overlap)
        );

        dV_table.resize(
            (axis_info_vr.get_num_grid()+2*overlap)*
            (axis_info_v_theta.get_num_grid()+2*overlap)*
            (axis_info_v_phi.get_num_grid()+2*overlap)
        );

        precompute_dV();
    }
    DistributionFunction(const DistributionFunction&r):
        axis_info_x(r.axis_info_x),
        axis_info_vr(r.axis_info_vr),
        axis_info_v_theta(r.axis_info_v_theta),
        axis_info_v_phi(r.axis_info_v_phi),
        overlap(overlap),
        f(r.f),
        r_theta_phi(r.r_theta_phi),
        theta_phi(r.theta_phi),
        phi(r.phi)
    {}
    inline const Value& at(
        Index id_x,
        Index id_vr,
        Index id_vtheta,
        Index id_vphi
    )const{
        return f
            [
                (id_x+overlap)        * r_theta_phi   +
                (id_vr+overlap)       * theta_phi     +
                (id_vtheta+overlap)   * phi           +
                (id_vphi+overlap)
            ];
    }
    inline Value& at(
        Index id_x,
        Index id_vr,
        Index id_vtheta,
        Index id_vphi
    ){
        return f
            [
                (id_x+overlap)        * r_theta_phi   +
                (id_vr+overlap)       * theta_phi     +
                (id_vtheta+overlap)   * phi           +
                (id_vphi+overlap)
            ];
    }

    void apply_boundary_condition(){
       
        //vr>vr_maxを０にする
        const Index x_size  = axis_info_x.get_num_grid();
        const Index r_size  = axis_info_vr.get_num_grid();
        const Index t_size  = axis_info_v_theta.get_num_grid();
        const Index p_size  = axis_info_v_phi.get_num_grid();

        // ============================
        // x: 周期境界
        // ============================
        // ループ順: r,t,p 外側、k は内側の短いループ
        #pragma omp parallel for collapse(3) schedule(static)
        for(Index r=0;r<r_size;++r){
            for(Index t=0;t<t_size;++t){
                for(Index p=0;p<p_size;++p){
                    // 左ゴーストセルに右端をコピー
                    for(Index k=1;k<=overlap;++k){
                        this->at(-k, r, t, p) = this->at(x_size-k, r, t, p);
                    }
                    // 右ゴーストセルに左端をコピー
                    for(Index k=0;k<overlap;++k){
                        this->at(x_size+k, r, t, p) = this->at(k, r, t, p);
                    }
                }
            }
        }

        // ============================
        // vr: 上端を 0
        // ============================
        // x,r,t,p 全体を大きくcollapseして並列化
        #pragma omp parallel for collapse(4) schedule(static)
        for(Index x=0;x<x_size;++x){
            for(Index r=r_size;r<r_size+overlap;++r){
                for(Index t=0;t<t_size;++t){
                    for(Index p=0;p<p_size;++p){
                        this->at(x,r,t,p)=0.0;
                    }
                }
            }
        }

        // ============================
        // v_r < 0  -> (-vr, pi-theta, phi+pi)
        // ============================
        #pragma omp parallel for collapse(4) schedule(static)
        for(Index x=0;x<x_size;++x){
            for(Index r=-overlap;r<0;++r){
                for(Index t=0;t<t_size;++t){
                    for(Index p=0;p<p_size;++p){
                        Index r_mapped = -r-1; // vr→-vr
                        Index t_mapped = t_size-1-t; // θ→π-θ
                        Index p_mapped = (p + p_size/2) % p_size; // φ→φ+π
                        this->at(x,r,t,p) = this->at(x,r_mapped,t_mapped,p_mapped);
                    }
                }
            }
        }

        // ============================
        // theta < 0 -> (vr, -theta, phi+PI)
        // ============================
        #pragma omp parallel for collapse(4) schedule(static)
        for(Index x=0;x<x_size;++x){
            for(Index r=0;r<r_size;++r){
                for(Index t=-overlap;t<0;++t){
                    for(Index p=0;p<p_size;++p){
                        Index t_mapped = -t-1;
                        Index p_mapped = (p + p_size/2) % p_size;
                        this->at(x,r,t,p) = this->at(x,r,t_mapped,p_mapped);
                    }
                }
            }
        }

        // ============================
        // theta > PI -> (vr, 2PI-theta, phi+PI)
        // ============================
        #pragma omp parallel for collapse(4) schedule(static)
        for(Index x=0;x<x_size;++x){
            for(Index r=0;r<r_size;++r){
                for(Index t=t_size;t<t_size+overlap;++t){
                    for(Index p=0;p<p_size;++p){
                        Index t_mapped = 2*t_size - t - 1;
                        Index p_mapped = (p + p_size/2) % p_size;
                        this->at(x,r,t,p) = this->at(x,r,t_mapped,p_mapped);
                    }
                }
            }
        }

        // ============================
        // phi < 0 -> phi + 2PI
        // ============================
        #pragma omp parallel for collapse(4) schedule(static)
        for(Index x=0;x<x_size;++x){
            for(Index r=0;r<r_size;++r){
                for(Index t=0;t<t_size;++t){
                    for(Index p=-overlap;p<0;++p){
                        Index p_mapped = p + p_size;
                        this->at(x,r,t,p) = this->at(x,r,t,p_mapped);
                    }
                }
            }
        }

        // ============================
        // phi > 2PI -> phi - 2PI
        // ============================
        #pragma omp parallel for collapse(4) schedule(static)
        for(Index x=0;x<x_size;++x){
            for(Index r=0;r<r_size;++r){
                for(Index t=0;t<t_size;++t){
                    for(Index p=p_size;p<p_size+overlap;++p){
                        Index p_mapped = p - p_size;
                        this->at(x,r,t,p) = this->at(x,r,t,p_mapped);
                    }
                }
            }
        }
    }

    void mul(Value val){
        Index f_size = f.size();
        for(Index i=0;i<f_size;++i){
            f[i]*=val;
        }
    }

    void add(const DistributionFunction& r_f){
        Index f_size = f.size();
        for(Index i=0;i<f_size;++i){
            f[i] += r_f.f[i];
        }
    }
};
#endif// DISTRIBUTION_FUNCTION_H