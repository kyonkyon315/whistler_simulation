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

using Value = double;

// calc_flux_rightward は「中心 f0 を用いる (fm2,fm1,f0,fp1,fp2)」で
// nyu >= 0 (左→右 情報) を想定した計算（先に示したものと同じ）
static inline Value calc_flux_rightward(
    Value fm2, Value fm1, Value f0, Value fp1, Value fp2,
    Value nyu
){
    Value f_max1 = std::max(std::max(fm1, f0),
                     std::min(2.0*fm1 - fm2, 2.0*f0 - fp1));
    Value f_max2 = std::max(std::max(fp1, f0),
                     std::min(2.0*fp1 - fp2, 2.0*f0 - fm1));
    Value f_min1 = std::min(std::min(fm1, f0),
                     std::max(2.0*fm1 - fm2, 2.0*f0 - fp1));
    Value f_min2 = std::min(std::min(fp1, f0),
                     std::max(2.0*fp1 - fp2, 2.0*f0 - fm1));

    Value f_max = std::max(f_max1, f_max2);
    Value f_min = std::min(f_min1, f_min2);

    Value L_pos;
    if (fp1 >= f0)
        L_pos = std::min(2.0*(f0 - f_min), (fp1 - f0));
    else
        L_pos = std::max(2.0*(f0 - f_max), (fp1 - f0));

    Value L_neg;
    if (f0 >= fm1)
        L_neg = std::min(2.0*(f_max - f0), (f0 - fm1));
    else
        L_neg = std::max(2.0*(f_min - f0), (f0 - fm1));

    Value U = (f0
            + (1.0 - nyu)*(2.0 - nyu)*L_pos
            + (1.0 - nyu)*(1.0 + nyu)*L_neg)*nyu;

    return U;
}

// 修正版：最悪ケースの7点を受け取り、符号に応じて正しいセンタリングを行う
Value finite_volume_method(
    Value f_im3, Value f_im2, Value f_im1,
    Value f_i,
    Value f_ip1, Value f_ip2, Value f_ip3,
    Value nyu
){
    // nyu = - v * dt / dx の定義をそのまま使う
    // delta_f = U_{i-1/2} - U_{i+1/2}
    Value delta_f_i = 0.0;

    if (nyu >= 0.0) {
        // 情報は左→右（upwind は 左 側。中心 f_i を使う）
        // U_{i+1/2} に対する中心 i と stencil i-2..i+2
        Value U_ip_half = calc_flux_rightward(
            /*fm2*/ f_im2,
            /*fm1*/ f_im1,
            /*f0*/  f_i,
            /*fp1*/ f_ip1,
            /*fp2*/ f_ip2,
            nyu
        );
        // U_{i-1/2} の中心は i-1 => stencil (i-3 .. i+1)
        Value U_im_half = calc_flux_rightward(
            /*fm2*/ f_im3,
            /*fm1*/ f_im2,
            /*f0*/  f_im1,
            /*fp1*/ f_i,
            /*fp2*/ f_ip1,
            nyu
        );
        delta_f_i = U_im_half - U_ip_half;
    } else {
        // nyu < 0: 情報は右→左（upwind は 右 側）。中心を i+1 にシフトして再構築する。
        Value abs_nyu = -nyu;

        // U_{i+1/2} の中心は i+1 -> stencil (i-1 .. i+3)
        Value U_ip_half = calc_flux_rightward(
            /*fm2*/ f_im1,  // (i+1)-2 = i-1
            /*fm1*/ f_i,    // (i+1)-1 = i
            /*f0*/  f_ip1,  // (i+1)
            /*fp1*/ f_ip2,  // (i+2)
            /*fp2*/ f_ip3,  // (i+3)
            abs_nyu
        );
        // U_{i-1/2} の中心は i -> stencil (i-2 .. i+2)
        Value U_im_half = calc_flux_rightward(
            /*fm2*/ f_im2,  // i-2
            /*fm1*/ f_im1,  // i-1
            /*f0*/  f_i,    // i
            /*fp1*/ f_ip1,  // i+1
            /*fp2*/ f_ip2,  // i+2
            abs_nyu
        );
        delta_f_i = U_im_half - U_ip_half;
    }

    return delta_f_i;
}

#endif //FINITE_VOLUME_METHOD_H