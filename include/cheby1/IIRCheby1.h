#pragma once

#include "common/IIRFilter.h"

namespace tiny_iir {

/**
 * @brief Chebyshev Type I filter
 *
 * @tparam N   Filter order.
 * @tparam T   Data type.
 * @tparam PASS_TYPE   Pass type (low, high).
 */
template<size_t N = 2, typename T = double, FilterPassType PASS_TYPE = FilterPassType::LOW_PASS>
class IIRCheby1 : public IIRFilter<N, T, PASS_TYPE> {
public:
    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::LOW_PASS
                                         || _PT == FilterPassType::HIGH_PASS)>>
    IIRCheby1(double normalized_cutoff_frequency, double passband_ripple_db);

    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::BAND_PASS
                                         || _PT == FilterPassType::BAND_STOP)>>
    IIRCheby1(double normalized_lowcut_freq, double normalized_highcut_freq, double passband_ripple_db);

    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::LOW_PASS
                                         || _PT == FilterPassType::HIGH_PASS)>>
    void configure(double normalized_cutoff_frequency, double passband_ripple_db);

    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::BAND_PASS
                                         || _PT == FilterPassType::BAND_STOP)>>
    void configure(double normalized_lowcut_freq, double normalized_highcut_freq, double passband_ripple_db);

private:
    void init(double ripple_db);

    PoleZeroPair get_pole_zero_pairs_s_plane(unsigned int i) final;

    PoleZeroPair get_pole_zero_real_axis() final;

    double _sinh_mu = 0;
    double _cosh_mu = 0;
};


template<size_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType _PT, typename>
IIRCheby1<N, T, PASS_TYPE>::IIRCheby1(double normalized_cutoff_frequency, double passband_ripple_db)
        : IIRFilter<N, T, PASS_TYPE>() {
    configure(normalized_cutoff_frequency, passband_ripple_db);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType _PT, typename>
IIRCheby1<N, T, PASS_TYPE>::IIRCheby1(double normalized_lowcut_freq, double normalized_highcut_freq,
                                      double passband_ripple_db)
        : IIRFilter<N, T, PASS_TYPE>() {
    configure(normalized_lowcut_freq, normalized_highcut_freq, passband_ripple_db);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType _PT, typename>
void IIRCheby1<N, T, PASS_TYPE>::configure(double normalized_cutoff_frequency, double passband_ripple_db) {
    init(passband_ripple_db);
    IIRFilter<N, T, PASS_TYPE>::calculate_cascades(normalized_cutoff_frequency);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType _PT, typename>
void
IIRCheby1<N, T, PASS_TYPE>::configure(double normalized_lowcut_freq, double normalized_highcut_freq,
                                      double passband_ripple_db) {
    init(passband_ripple_db);
    IIRFilter<N, T, PASS_TYPE>::calculate_cascades(normalized_lowcut_freq, normalized_highcut_freq);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
PoleZeroPair IIRCheby1<N, T, PASS_TYPE>::get_pole_zero_pairs_s_plane(unsigned int i) {
    const double phi = (2 * i + 1) * M_PI_2 / N; // Angle from the imaginary axis
    const double sin_phi = std::sin(phi);
    const double cos_phi = std::cos(phi);

    const double pole_s_real = -_sinh_mu * sin_phi;
    const double pole_s_imag = _cosh_mu * cos_phi;
    return {Complex{pole_s_real, pole_s_imag},
            Complex{IIRFilter<N, T, PASS_TYPE>::INFINITY_VALUE, 0}};
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
PoleZeroPair IIRCheby1<N, T, PASS_TYPE>::get_pole_zero_real_axis() {
    return {Complex{-_sinh_mu, 0},
            Complex{IIRFilter<N, T, PASS_TYPE>::INFINITY_VALUE, 0}};
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
void IIRCheby1<N, T, PASS_TYPE>::init(double ripple_db) {
    ripple_db = std::abs(ripple_db);

    if constexpr (N & 1) {
        IIRFilter<N, T, PASS_TYPE>::_gain_double = 1.0;
    } else {
        IIRFilter<N, T, PASS_TYPE>::_gain_double = std::exp(-ripple_db / 20 * M_LN10);
    }

    const double epsilon = std::sqrt(std::exp(ripple_db * 0.1 * M_LN10) - 1);
    const double mu = std::asinh(1.0 / epsilon) / N;
    _sinh_mu = std::sinh(mu);
    _cosh_mu = std::cosh(mu);
}

}