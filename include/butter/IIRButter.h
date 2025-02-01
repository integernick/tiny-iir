#pragma once

#include "common/IIRFilter.h"

namespace tiny_iir {
template<size_t N = 2, typename T = double, FilterPassType PASS_TYPE = FilterPassType::LOW_PASS>
class IIRButter : public IIRFilter<N, T, PASS_TYPE> {
public:
    IIRButter(double normalized_cutoff_frequency);

    PoleZeroPair get_pole_zero_pairs_s_plane(unsigned int i) final;

    PoleZeroPair get_pole_zero_real_axis() final;

    void configure(double Wn);

private:
    double _d_phi = 0;
};

template<size_t N, typename T, FilterPassType PASS_TYPE>
IIRButter<N, T, PASS_TYPE>::IIRButter(double normalized_cutoff_frequency) {
    configure(normalized_cutoff_frequency);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
PoleZeroPair IIRButter<N, T, PASS_TYPE>::get_pole_zero_pairs_s_plane(unsigned int i) {
    const double phi = (2 * i + 1) * _d_phi;
    const double pole_s_real = -std::sin(phi);
    const double pole_s_imag = std::cos(phi);
    return {{pole_s_real,                                pole_s_imag},
            {IIRFilter<N, T, PASS_TYPE>::INFINITY_VALUE, 0}};
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
PoleZeroPair IIRButter<N, T, PASS_TYPE>::get_pole_zero_real_axis() {
    return {-1.0,
            IIRFilter<N, T, PASS_TYPE>::INFINITY_VALUE};
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
void IIRButter<N, T, PASS_TYPE>::configure(double Wn) {
    IIRFilter<N, T, PASS_TYPE>::_gain_double = 1.0;
    _d_phi = M_PI_2 / N;
    IIRFilter<N, T, PASS_TYPE>::configure(Wn);
}

}