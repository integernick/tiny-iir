#pragma once

#include "common/IIRFilter.h"

namespace tiny_iir {

template<uint32_t N = 2, 
        typename T = double, 
        FilterPassType PASS_TYPE = FilterPassType::LOW_PASS,
        uint32_t CROSSFADE_SAMPLES = 0>
class IIRButter : public IIRFilter<N, T, PASS_TYPE, CROSSFADE_SAMPLES> {
public:
    using IIRFilterBase = IIRFilter<N, T, PASS_TYPE, CROSSFADE_SAMPLES>;

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
        init_analog();
        IIRFilterBase::calculate_cascades(normalized_cutoff_frequency);
    }

    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::BAND_PASS
                                         || _PT == FilterPassType::BAND_STOP)>>
    void configure(double normalized_lowcut_freq, double normalized_highcut_freq) {
        init_analog();
        IIRFilterBase::calculate_cascades(normalized_lowcut_freq, normalized_highcut_freq);
    }

private:
    void init_analog();
};

template<uint32_t N, typename T, FilterPassType PASS_TYPE, uint32_t CROSSFADE_SAMPLES>
template<FilterPassType _PT, typename>
IIRButter<N, T, PASS_TYPE, CROSSFADE_SAMPLES>::IIRButter(double normalized_lowcut_freq, double normalized_highcut_freq)
        : IIRFilterBase() {
    configure(normalized_lowcut_freq, normalized_highcut_freq);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, uint32_t CROSSFADE_SAMPLES>
template<FilterPassType _PT, typename>
IIRButter<N, T, PASS_TYPE, CROSSFADE_SAMPLES>::IIRButter(double normalized_cutoff_frequency)
        : IIRFilterBase() {
    configure(normalized_cutoff_frequency);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, uint32_t CROSSFADE_SAMPLES>
void IIRButter<N, T, PASS_TYPE, CROSSFADE_SAMPLES>::init_analog() {
    IIRFilterBase::_gain_double = 1.0;

    if constexpr (N & 1) {
        IIRFilterBase::_analog_pole_zero_pairs[(N + 1) / 2 - 1] = {
                -1.0,
                IIRFilterBase::INFINITY_VALUE
        };
    }

    for (uint32_t i = 0; i < N / 2; ++i) {
        const double phi = static_cast<double>(2 * i + 1) * M_PI_2 / N; // Angle from the imaginary axis
        const double pole_s_real = -std::sin(phi);
        const double pole_s_imag = std::cos(phi);
        IIRFilterBase::_analog_pole_zero_pairs[i] = {
                {pole_s_real, pole_s_imag},
                {INFINITY_VALUE, 0}
        };
    }
}

}