#ifndef FINITE_VOLUME_METHOD_H
#define FINITE_VOLUME_METHOD_H
using Value = double;

/*

[Takayuki Umeda et.al 2008, "A conservative and non-oscillatory scheme for Vlasov code simulations"]

df/dt + v df/dx = 0
を解くルーチンをここに書く
nyu = - v Δt/Δx
とすると、

f(i,t+Δt) = f(i,t) + U_(i-1/2)(nyu) - U_(i+1/2)(nyu)

ただし、
U_(i+1/2)(nyu) = nyu f(i) + nyu(1-nyu)(2-nyu)L_i_positive
                            + nyu(1-nyu)(1+nyu)L_i_negative

L_i_positive =  min[2(f_i - f_min), (f_i+1 - f_i)] if f_i+1 >= f_i
                max[2(f_i - f_max), (f_i+1 - f_i)] else

L_i_negative =  min[2(f_max - f_i), (f_i - f_i-1)] if f_i >= f_i-1
                max[2(f_min - f_i), (f_i - f_i-1)] else

f_max = max[f_max1, f_max2]
f_min = min[f_min1, f_min2]

f_max1 = max[max[f_i-1, f_i], min[2f_i-1 - f_i-2, 2f_i - f_i+1]]
f_max2 = max[max[f_i+1, f_i], min[2f_i+1 - f_i+2, 2f_i - f_i-1]]
f_min1 = min[min[f_i-1, f_i], max[2f_i-1 - f_i-2, 2f_i - f_i+1]]
f_min2 = min[min[f_i+1, f_i], max[2f_i+1 - f_i+2, 2f_i - f_i-1]]

*/

#include <algorithm> // for std::min/std::max

Value finite_volume_method(
    Value f_i_minus_3,
    Value f_i_minus_2,
    Value f_i_minus_1,
    Value f_i,
    Value f_i_plus_1,
    Value f_i_plus_2,
    Value nyu
){
    auto U_at_half = [&](Value fm2, Value fm1, Value f0, Value fp1, Value fp2) -> Value {
        // fm2 = f_{i-2}, fm1 = f_{i-1}, f0 = f_i, fp1 = f_{i+1}, fp2 = f_{i+2}
        // compute f_max1, f_max2, f_min1, f_min2
        Value f_max1 = std::max(std::max(fm1, f0), std::min(2*fm1 - fm2, 2*f0 - fp1));
        Value f_max2 = std::max(std::max(fp1, f0), std::min(2*fp1 - fp2, 2*f0 - fm1));
        Value f_min1 = std::min(std::min(fm1, f0), std::max(2*fm1 - fm2, 2*f0 - fp1));
        Value f_min2 = std::min(std::min(fp1, f0), std::max(2*fp1 - fp2, 2*f0 - fm1));

        Value f_max = std::max(f_max1, f_max2);
        Value f_min = std::min(f_min1, f_min2);

        // L_positive
        Value L_pos;
        if (fp1 >= f0) {
            L_pos = std::min( 2*(f0 - f_min), (fp1 - f0) );
        } else {
            L_pos = std::max( 2*(f0 - f_max), (fp1 - f0) );
        }

        // L_negative
        Value L_neg;
        if (f0 >= fm1) {
            L_neg = std::min( 2*(f_max - f0), (f0 - fm1) );
        } else {
            L_neg = std::max( 2*(f_min - f0), (f0 - fm1) );
        }

        // flux U_{i+1/2}(nyu) (or general center-based half flux)
        Value U = nyu * f0
                + nyu*(1 - nyu)*(2 - nyu) * L_pos
                + nyu*(1 - nyu)*(1 + nyu) * L_neg;
        return U;
    };

    // U_{i+1/2} uses center f_i and neighbors f_{i-2}..f_{i+2}
    Value U_ip_half = U_at_half(f_i_minus_2, f_i_minus_1, f_i, f_i_plus_1, f_i_plus_2);

    // U_{i-1/2} uses center f_{i-1} and neighbors f_{i-3}..f_{i+1}
    Value U_im_half = U_at_half(f_i_minus_3, f_i_minus_2, f_i_minus_1, f_i, f_i_plus_1);

    // delta_f = U_{i-1/2} - U_{i+1/2}
    Value delta_f_i = U_im_half - U_ip_half;

    return delta_f_i;
}


#endif //FINITE_VOLUME_METHOD_H