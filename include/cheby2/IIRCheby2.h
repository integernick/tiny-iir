#pragma once

#include "common/IIRFilter.h"

namespace tiny_iir {
template<size_t N = 2, typename T = double, FilterPassType PASS_TYPE = FilterPassType::LOW_PASS>
class IIRCheby2 : public IIRFilter<N, T, PASS_TYPE> {
public:
    IIRCheby2(double normalized_cutoff_frequency, double ripple_db);

    PoleZeroPair get_pole_zero_pairs_s_plane(unsigned int i) final;

    PoleZeroPair get_pole_zero_real_axis() final;

    void configure(double Wn, double ripple_db);

private:
    double _d_phi = 0;
    double _sinh_mu = 0;
    double _cosh_mu = 0;
};


template<size_t N, typename T, FilterPassType PASS_TYPE>
IIRCheby2<N, T, PASS_TYPE>::IIRCheby2(double normalized_cutoff_frequency, double ripple_db) {
    configure(normalized_cutoff_frequency, ripple_db);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
PoleZeroPair IIRCheby2<N, T, PASS_TYPE>::get_pole_zero_pairs_s_plane(unsigned int i) {
    const double phi = (2 * i + 1) * _d_phi;
    const double sin_phi = std::sin(phi);
    const double cos_phi = std::cos(phi);

    const double pole_s_real = -_sinh_mu * sin_phi;
    const double pole_s_imag = _cosh_mu * cos_phi;
    const double pole_squared = pole_s_real * pole_s_real + pole_s_imag * pole_s_imag;
    const double zero_imag = 1 / cos_phi;
    return {std::complex<double>{pole_s_real / pole_squared, pole_s_imag / pole_squared},
            std::complex<double>{0, zero_imag}};
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
PoleZeroPair IIRCheby2<N, T, PASS_TYPE>::get_pole_zero_real_axis() {
    return {-1.0 / _sinh_mu,
            IIRFilter<N, T, PASS_TYPE>::INFINITY_VALUE};
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
void IIRCheby2<N, T, PASS_TYPE>::configure(double Wn, double ripple_db) {
    IIRFilter<N, T, PASS_TYPE>::_cascade_filter.set_gain(1.0);

    const double delta = 1.0 / std::sqrt(std::exp(ripple_db * 0.1 * M_LN10) - 1);
    const double mu = std::asinh(1.0 / delta) / N;
    _d_phi = M_PI_2 / N;
    _sinh_mu = std::sinh(mu);
    _cosh_mu = std::cosh(mu);

    IIRFilter<N, T, PASS_TYPE>::configure(Wn);
}

}