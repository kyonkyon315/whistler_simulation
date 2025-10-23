#ifndef AXISINFO_H
#define AXISINFO_H
#include <cassert>
#define PI 3.14159265358979323846
class AxisInfo_x{
    using Value = double;
    using Index = int;
    private:
    const Value begin;
    const Value end;
    const Value grid_size;
    const Index num_grid;
    public:
    AxisInfo_x(Value begin,Value end,Index num_grid):
        begin(begin),
        end(end),
        grid_size((end-begin)/static_cast<Value>(num_grid-1)),
        num_grid(num_grid)
    {}
    Value operator[](Index id){return begin+grid_size*static_cast<Value>(id);}
    Value get_begin()const{return begin;}
    Value get_end()const{return end;}
    Index get_num_grid()const{return num_grid;}
};

class AxisInfo_r{
    using Value = double;
    using Index = int;
    private:
    const Value begin;
    const Value end;
    const Value grid_size;
    const Index num_grid;

    public:

    //0除算を避けるためにbegin=grid_size/2とする
    AxisInfo_r(Value end,Index num_grid):
        begin(end/static_cast<Value>(2*num_grid-1)),
        end(end),
        grid_size(2.*begin),
        num_grid(num_grid)
    {}
    Value operator[](Index id){return begin+grid_size*static_cast<Value>(id);}
    Value get_begin()const{return begin;}
    Value get_end()const{return end;}
    Index get_num_grid()const{return num_grid;}
};

class AxisInfo_theta{
    using Value = double;
    using Index = int;
    private:
    const Value begin;
    const Value end;
    const Value grid_size;
    const Index num_grid;

    public:

    //0除算を避けるためにbegin=grid_size/2とする
    AxisInfo_theta(Index num_grid):
        begin(PI/static_cast<Value>(2*num_grid)),
        end(PI-begin),
        grid_size(2.*begin),
        num_grid(num_grid)
    {}
    Value operator[](Index id){return begin+grid_size*static_cast<Value>(id);}
    Value get_begin()const{return begin;}
    Value get_end()const{return end;}
    Index get_num_grid()const{return num_grid;}
};

class AxisInfo_phi{
    using Value = double;
    using Index = int;
    private:
    const Value begin;
    const Value end;
    const Value grid_size;
    const Index num_grid;

    public:

    //0除算を避けるためにbegin=grid_size/2とする
    AxisInfo_phi(Index num_grid):
        begin(0.),
        grid_size(2.*PI/static_cast<Value>(num_grid)),
        end(2.*PI-grid_size),
        num_grid(num_grid)
    {
        assert(num_grid%2==0);
    }
    Value operator[](Index id){return begin+grid_size*static_cast<Value>(id);}
    Value get_begin()const{return begin;}
    Value get_end()const{return end;}
    Index get_num_grid()const{return num_grid;}
};

#endif// AXISINFO_H
