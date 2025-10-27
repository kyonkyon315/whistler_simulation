#ifndef DISTRIBUTION_FUNCTION_DEVELOPPER_H
#define DISTRIBUTION_FUNCTION_DEVELOPPER_H
#include "distribution_function.h"
#include "current.h"
#include "electro_magnetic_field.h"
#include "finite_volume_method.h"
#include <string>

using Value = double;
using Index = int;

class DistributionFunctionDevelopper{
private:
    const AxisInfo_x& axis_x;
    const AxisInfo_r& axis_r;
    const AxisInfo_theta& axis_t;
    const AxisInfo_phi& axis_p;
public:

    DistributionFunctionDevelopper(
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

    //fをv・∇fの寄与に基づきdt/2発展
    void advect_x(DistributionFunction &f,Value dt){
        
        // まず境界（ゴーストセル）を最新化しておく
        f.apply_boundary_condition();

        const Index x_size   = axis_x.get_num_grid();
        const Index r_size   = axis_r.get_num_grid();
        const Index t_size   = axis_t.get_num_grid();
        const Index p_size   = axis_p.get_num_grid();
        const Value dx       = axis_x.get_grid_size();


        // 速度セルごとの x 軸一時バッファ
        std::vector<Value> new_x(x_size);

        // 並列化: 速度空間ループを外側にし、x は短いループ（書き戻しは各スレッド別領域）
        #pragma omp parallel for collapse(3) schedule(static)
        for(Index r=0;r<r_size;++r){
            for(Index t=0;t<t_size;++t){
                for(Index p=0;p<p_size;++p){

                    // 速度（x成分）を取得（この速度はこの (r,t,p) で定数）
                    Value v_x = f.get_velocity_at<Direction::x>(0, r, t, p);

                    // Courant 数 nyu = - v * dt / dx（元の式に合わせる）
                    Value nyu = - v_x * dt / dx;

                    // 各 x 点に対して有限体積差分（7点ステンシル）
                    for(Index ix = 0; ix < x_size; ++ix){
                        // stencil indices: i-3..i+2 relative to interior indexing [0..x_size-1]
                        // 注意: DistributionFunction::at はゴーストセルのアクセスも受け付ける（overlap あり）
                        Value f_i_minus_3 = f.at(ix - 3, r, t, p);
                        Value f_i_minus_2 = f.at(ix - 2, r, t, p);
                        Value f_i_minus_1 = f.at(ix - 1, r, t, p);
                        Value f_i         = f.at(ix    , r, t, p);
                        Value f_i_plus_1  = f.at(ix + 1, r, t, p);
                        Value f_i_plus_2  = f.at(ix + 2, r, t, p);
                        Value f_i_plus_3  = f.at(ix + 3, r, t, p);

                        Value delta = finite_volume_method(
                            f_i_minus_3, f_i_minus_2, f_i_minus_1,
                            f_i, 
                            f_i_plus_1, f_i_plus_2, f_i_plus_3,
                            nyu
                        );

                        new_x[ix] = f_i + delta;
                    }

                    // 書き戻し（内側点のみ）
                    for(Index ix = 0; ix < x_size; ++ix){
                        f.at(ix, r, t, p) = new_x[ix];
                    }
                }
            }
        }
    }

    //fをqE・∇_v fの寄与に基づきdt/2発展
    void advect_E(DistributionFunction &f,const ElectroMagneticField &,Value dt){
    }

    //fをqVxB・∇_v fの寄与に基づきdt発展
    void advect_B(DistributionFunction &f,const ElectroMagneticField &,Value dt){

    }

};

#endif// DISTRIBUTION_FUNCTION_DEVELOPPER_H
