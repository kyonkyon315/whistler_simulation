#ifndef CURRENT_H
#define CURRENT_H
#include "axis_info.h"
#include "direction.h"

#include <vector>


class Current{
    using Value = double;
    using Index = int;
    private:
    const AxisInfo_x& axis_info_x;
    const Index grid_size;
    const Index overlap;
    const Index stride;
    std::vector<Value> current;
    
    public:

    //例えばoverlap=5の場合、.at(-5)や.at(grid_size+4)が許可される。.at(-6)は未定義動作
    Current(const AxisInfo_x& axis_info_x,Index overlap=5):
        axis_info_x(axis_info_x),
        overlap(overlap),
        grid_size(axis_info_x.get_num_grid()),
        stride(grid_size+overlap*2)
    {
        current.resize(3*(stride));
    }

    template<Direction dir>
    const Value& at(Index id) const {
        constexpr Index d       = (dir==Direction::x ? 0 : (dir==Direction::y ? 1 : 2));
        constexpr Index idx     = d ;
        return current[idx*stride + id + overlap];
    }
    template<Direction dir>
    Value& at(Index id){
        constexpr Index d       = (dir==Direction::x ? 0 : (dir==Direction::y ? 1 : 2));
        constexpr Index idx     = d ;
        return current[idx*stride + id + overlap];
    }


    void apply_boundary_condition() {
    // 全ての成分（Ex,Ey,Ez,Bx,By,Bz の main/half）について処理
    for (int comp = 0; comp < 12; ++comp) {

        // 左側のゴーストセルに右端の値をコピー
        for (Index i = 0; i < overlap; ++i) {
            at<Direction::x>(grid_size+i) = at<Direction::x>(i);  
            at<Direction::y>(grid_size+i) = at<Direction::y>(i);  
            at<Direction::z>(grid_size+i) = at<Direction::z>(i);  
        }

        // 右側のゴーストセルに左端の値をコピー
        for (Index i = 1; i <= overlap; ++i) {
            at<Direction::x>(-i) = at<Direction::x>(grid_size-i);  
            at<Direction::y>(-i) = at<Direction::y>(grid_size-i);  
            at<Direction::z>(-i) = at<Direction::z>(grid_size-i);  
        }
    }
}
    
};
#endif// CURRENT_H