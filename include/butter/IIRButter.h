#pragma once

#include "common/IIRFilter.h"

namespace tiny_iir {
template<size_t N = 2, typename T = double, FilterPassType PASS_TYPE = FilterPassType::LOW_PASS>
class IIRButter : public IIRFilter<N, T, PASS_TYPE> {
public:
    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::LOW_PASS
                                         || _PT == FilterPassType::HIGH_PASS)>>
    IIRButter(double normalized_cutoff_frequency);

    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::BAND_PASS
                                         || _PT == FilterPassType::BAND_STOP)>>
    IIRButter(double normalized_lowcut_freq, double normalized_highcut_freq);

    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::LOW_PASS
                                         || _PT == FilterPassType::HIGH_PASS)>>
    void configure(double normalized_cutoff_frequency) {
        init();
        IIRFilter<N, T, PASS_TYPE>::calculate_cascades(normalized_cutoff_frequency);
    }

    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::BAND_PASS
                                         || _PT == FilterPassType::BAND_STOP)>>
    void configure(double normalized_lowcut_freq, double normalized_highcut_freq) {
        init();
        IIRFilter<N, T, PASS_TYPE>::calculate_cascades(normalized_lowcut_freq, normalized_highcut_freq);
    }

private:
    void init();

    PoleZeroPair get_pole_zero_pairs_s_plane(unsigned int i) final;

    PoleZeroPair get_pole_zero_real_axis() final;

    double _d_phi = 0;
};

template<size_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType _PT, typename>
IIRButter<N, T, PASS_TYPE>::IIRButter(double normalized_lowcut_freq, double normalized_highcut_freq)
        : IIRFilter<N, T, PASS_TYPE>() {
    configure(normalized_lowcut_freq, normalized_highcut_freq);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType _PT, typename>
IIRButter<N, T, PASS_TYPE>::IIRButter(double normalized_cutoff_frequency)
        : IIRFilter<N, T, PASS_TYPE>() {
    configure(normalized_cutoff_frequency);
}


template<size_t N, typename T, FilterPassType PASS_TYPE>
PoleZeroPair IIRButter<N, T, PASS_TYPE>::get_pole_zero_pairs_s_plane(unsigned int i) {
    const double phi = (2 * i + 1) * _d_phi; // Angle from the imaginary axis
    const double pole_s_real = -std::sin(phi);
    const double pole_s_imag = std::cos(phi);
    return {{pole_s_real, pole_s_imag},
            {IIRFilter<N, T, PASS_TYPE>::INFINITY_VALUE, 0}};
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
PoleZeroPair IIRButter<N, T, PASS_TYPE>::get_pole_zero_real_axis() {
    return {-1.0,
            IIRFilter<N, T, PASS_TYPE>::INFINITY_VALUE};
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
void IIRButter<N, T, PASS_TYPE>::init() {
    IIRFilter<N, T, PASS_TYPE>::_gain_double = 1.0;
    _d_phi = M_PI_2 / N;
}

}