#ifndef ELECTROMAGNETICFIELD_H
#define ELECTROMAGNETICFIELD_H
#include "axis_info.h"
#include <vector>

enum class FieldType { E, B };
enum class Direction { x, y, z };
enum class TimeStagger { main, half };

class ElectroMagneticField{
    using Value = double;
    using Index = int;
    private:
    const AxisInfo_x& axis_info_x;
    const Index grid_size;
    const Index overlap;
    const Index stride;
    std::vector<Value> field;
    
    public:

    //例えばoverlap=5の場合、.at(-5)や.at(grid_size+4)が許可される。.at(-6)は未定義動作
    ElectroMagneticField(const AxisInfo_x& axis_info_x,Index overlap=5):
        axis_info_x(axis_info_x),
        overlap(overlap),
        grid_size(axis_info_x.get_num_grid()),
        stride(grid_size+overlap*2)
    {
        field.resize(12*(stride));
    }

    template<FieldType e_or_b, Direction dir, TimeStagger main_or_half>
    const Value& at(Index id) const {
        constexpr Index base    = (e_or_b==FieldType::E ? 0 : 6);
        constexpr Index d       = (dir==Direction::x ? 0 : (dir==Direction::y ? 1 : 2));
        constexpr Index m_or_h  =(main_or_half==TimeStagger::main ? 0 : 3);
        constexpr Index idx     = base + d + m_or_h;
        return field[idx*stride + id + overlap];
    }
    template<FieldType e_or_b, Direction dir, TimeStagger main_or_half>
    Value& at(Index id){
        constexpr Index base    = (e_or_b==FieldType::E ? 0 : 6);
        constexpr Index d       = (dir==Direction::x ? 0 : (dir==Direction::y ? 1 : 2));
        constexpr Index m_or_h  =(main_or_half==TimeStagger::main ? 0 : 3);
        constexpr Index idx     = base + d + m_or_h;
        return field[idx*stride + id + overlap];
    }

    template<FieldType e_or_b, Direction dir>
    const Value& at(Index id) const {
        constexpr Index base    = (e_or_b==FieldType::E ? 0 : 6);
        constexpr Index d       = (dir==Direction::x ? 0 : (dir==Direction::y ? 1 : 2));
        constexpr Index m_or_h  = 0;
        constexpr Index idx     = base + d + m_or_h;
        return field[idx*stride + id + overlap];
    }
    template<FieldType e_or_b, Direction dir>
    Value& at(Index id){
        constexpr Index base    = (e_or_b==FieldType::E ? 0 : 6);
        constexpr Index d       = (dir==Direction::x ? 0 : (dir==Direction::y ? 1 : 2));
        constexpr Index m_or_h  = 0;
        constexpr Index idx     = base + d + m_or_h;
        return field[idx*stride + id + overlap];
    }

    void apply_boundary_condition() {
    // 全ての成分（Ex,Ey,Ez,Bx,By,Bz の main/half）について処理
    for (int comp = 0; comp < 12; ++comp) {
        Value* f = field.data() + comp * stride;

        // 左側のゴーストセルに右端の値をコピー
        for (Index i = 0; i < overlap; ++i) {
            f[i] = f[grid_size + i];  
        }

        // 右側のゴーストセルに左端の値をコピー
        for (Index i = 0; i < overlap; ++i) {
            f[grid_size + overlap + i] = f[overlap + i];
        }
    }
}
    
};
#endif// ELECTROMAGNETICFIELD_H