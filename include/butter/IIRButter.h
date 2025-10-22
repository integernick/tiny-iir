#pragma once

#include "common/IIRFilter.h"

namespace tiny_iir {

template<uint32_t N = 2,
        typename T = double,
        FilterPassType PASS_TYPE = FilterPassType::LOW_PASS>
class IIRButter final : public IIRFilter<N, T, PASS_TYPE> {
public:
    using IIRFilterBase = IIRFilter<N, T, PASS_TYPE>;

    /**
     * @brief   Constructor (low-pass and high-pass).
     *
     * @param crossfade_samples  The number of samples to smooth the transition between old and new coefficients.
     */
    explicit IIRButter(double normalized_cutoff_frequency, uint32_t crossfade_samples = 0) requires (
    PASS_TYPE == FilterPassType::LOW_PASS or PASS_TYPE == FilterPassType::HIGH_PASS);

    /**
     * @brief   Constructor (band-pass and band-stop).
     *
     * @param crossfade_samples  The number of samples to smooth the transition between old and new coefficients.
     */
    IIRButter(double normalized_lowcut_freq, double normalized_highcut_freq, uint32_t crossfade_samples = 0) requires (
    PASS_TYPE == FilterPassType::BAND_PASS or PASS_TYPE == FilterPassType::BAND_STOP);

    /**
     * @brief   Configure the filter (low-pass and high-pass).
     *
     * @param normalized_cutoff_frequency  Normalized cutoff frequency.
     */
    void configure(double normalized_cutoff_frequency) requires (
    PASS_TYPE == FilterPassType::LOW_PASS or PASS_TYPE == FilterPassType::HIGH_PASS) {
        init_analog();
        IIRFilterBase::calculate_cascades(normalized_cutoff_frequency);
    }

    /**
     * @brief   Configure the filter (band-pass and band-stop).
     *
     * @param normalized_lowcut_freq  Normalized low-pass cutoff frequency.
     * @param normalized_highcut_freq  Normalized high-pass cutoff frequency.
     */
    void configure(double normalized_lowcut_freq, double normalized_highcut_freq) requires (
    PASS_TYPE == FilterPassType::BAND_PASS or PASS_TYPE == FilterPassType::BAND_STOP) {
        init_analog();
        IIRFilterBase::calculate_cascades(normalized_lowcut_freq, normalized_highcut_freq);
    }

private:
    /**
     * @brief   Get the analog gain.
     *
     * @return  The analog gain.
     */
    [[nodiscard]] double get_analog_gain() const final;

    /**
     * @brief   Initialize analog filter poles and zeros.
     */
    void init_analog() final;
};


template<uint32_t N, typename T, FilterPassType PASS_TYPE>
IIRButter<N, T, PASS_TYPE>::IIRButter(double normalized_cutoff_frequency, uint32_t crossfade_samples) requires (
PASS_TYPE == FilterPassType::LOW_PASS or PASS_TYPE == FilterPassType::HIGH_PASS)
        : IIRFilterBase(crossfade_samples) {
    configure(normalized_cutoff_frequency);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
IIRButter<N, T, PASS_TYPE>::IIRButter(double normalized_lowcut_freq, double normalized_highcut_freq,
                                      uint32_t crossfade_samples) requires (
PASS_TYPE == FilterPassType::BAND_PASS or PASS_TYPE == FilterPassType::BAND_STOP)
        : IIRFilterBase(crossfade_samples) {
    configure(normalized_lowcut_freq, normalized_highcut_freq);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
double IIRButter<N, T, PASS_TYPE>::get_analog_gain() const {
    return 1.0;
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
void IIRButter<N, T, PASS_TYPE>::init_analog() {
    if constexpr (N & 1) {
        IIRFilterBase::_analog_pole_zero_pairs[(N + 1) / 2 - 1] = {
                -1.0,
                INFINITY_VALUE
        };
    }

    for (uint32_t i = 0; i < N / 2; ++i) {
        const double phi = static_cast<double>(2 * i + 1) * M_PI_2 / N; // Angle from the imaginary axis
        const double pole_s_real = -std::sin(phi);
        const double pole_s_imag = std::cos(phi);
        IIRFilterBase::_analog_pole_zero_pairs[i] = {
                {pole_s_real,    pole_s_imag},
                {INFINITY_VALUE, 0}
        };
    }
}

}