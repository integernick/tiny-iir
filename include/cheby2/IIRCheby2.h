#pragma once

#include "common/IIRFilter.h"

namespace tiny_iir {
template<size_t N = 2, typename T = double, FilterPassType PASS_TYPE = FilterPassType::LOW_PASS>
class IIRCheby2 : public IIRFilter<N, T, PASS_TYPE> {
public:
    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::LOW_PASS
                                         || _PT == FilterPassType::HIGH_PASS)>>
    IIRCheby2(double normalized_cutoff_frequency, double stopband_ripple_db);

    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::BAND_PASS
                                         || _PT == FilterPassType::BAND_STOP)>>
    IIRCheby2(double normalized_lowcut_freq, double normalized_highcut_freq, double stopband_ripple_db);

    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::LOW_PASS
                                         || _PT == FilterPassType::HIGH_PASS)>>
    void configure(double normalized_cutoff_frequency, double stopband_ripple_db);

    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::BAND_PASS
                                         || _PT == FilterPassType::BAND_STOP)>>
    void configure(double normalized_lowcut_freq, double normalized_highcut_freq, double stopband_ripple_db);

private:
    void init_analog() final;

    double _stopband_ripple_db = 0;
};


template<size_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType _PT, typename>
IIRCheby2<N, T, PASS_TYPE>::IIRCheby2(double normalized_cutoff_frequency, double stopband_ripple_db)
        : IIRFilter<N, T, PASS_TYPE>() {
    configure(normalized_cutoff_frequency, stopband_ripple_db);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType _PT, typename>
IIRCheby2<N, T, PASS_TYPE>::IIRCheby2(double normalized_lowcut_freq, double normalized_highcut_freq,
                                      double stopband_ripple_db)
        : IIRFilter<N, T, PASS_TYPE>() {
    configure(normalized_lowcut_freq, normalized_highcut_freq, stopband_ripple_db);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType _PT, typename>
void IIRCheby2<N, T, PASS_TYPE>::configure(double normalized_cutoff_frequency, double stopband_ripple_db) {
    stopband_ripple_db = std::abs(stopband_ripple_db);
    if (stopband_ripple_db != _stopband_ripple_db) {
        _stopband_ripple_db = stopband_ripple_db;
        init_analog();
    }
    IIRFilter<N, T, PASS_TYPE>::calculate_cascades(normalized_cutoff_frequency);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType _PT, typename>
void
IIRCheby2<N, T, PASS_TYPE>::configure(double normalized_lowcut_freq, double normalized_highcut_freq,
                                      double stopband_ripple_db) {
    stopband_ripple_db = std::abs(stopband_ripple_db);
    if (stopband_ripple_db != _stopband_ripple_db) {
        _stopband_ripple_db = stopband_ripple_db;
        init_analog();
    }
    IIRFilter<N, T, PASS_TYPE>::calculate_cascades(normalized_lowcut_freq, normalized_highcut_freq);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
void IIRCheby2<N, T, PASS_TYPE>::init_analog() {
    IIRFilter<N, T, PASS_TYPE>::_cascade_filter.set_gain(1.0);

    const double delta = 1.0 / std::sqrt(std::exp(_stopband_ripple_db * 0.1 * M_LN10) - 1);
    const double mu = std::asinh(1.0 / delta) / N;
    const double sinh_mu = std::sinh(mu);
    const double cosh_mu = std::cosh(mu);

    if constexpr (N & 1) {
        IIRFilter<N, T, PASS_TYPE>::_analog_pole_zero_pairs[(N + 1) / 2 - 1] = {
                -1.0 / sinh_mu,
                IIRFilter<N, T, PASS_TYPE>::INFINITY_VALUE
        };
    }

    for (size_t i = 0; i < N / 2; ++i) {
        const double phi = static_cast<double>(2 * i + 1) * M_PI_2 / N; // Angle from the imaginary axis
        const double sin_phi = std::sin(phi);
        const double cos_phi = std::cos(phi);

        const double pole_s_real = -sinh_mu * sin_phi;
        const double pole_s_imag = cosh_mu * cos_phi;
        const double pole_squared = pole_s_real * pole_s_real + pole_s_imag * pole_s_imag;
        const double zero_imag = 1.0 / cos_phi;
        IIRFilter<N, T, PASS_TYPE>::_analog_pole_zero_pairs[i] = {
                Complex{pole_s_real / pole_squared, pole_s_imag / pole_squared},
                Complex{0, zero_imag}
        };
    }
}

}