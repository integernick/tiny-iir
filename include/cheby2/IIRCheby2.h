#pragma once

#include "common/IIRFilter.h"

namespace tiny_iir {
template<uint32_t N = 2, typename T = double, FilterPassType PASS_TYPE = FilterPassType::LOW_PASS>
class IIRCheby2 final : public IIRFilter<N, T, PASS_TYPE> {
public:
    /**
     * @brief   Constructor (low-pass and high-pass).
     *
     * @param crossfade_samples  The number of samples to smooth the transition between old and new coefficients.
     */
    IIRCheby2(double normalized_cutoff_frequency, double stopband_ripple_db, uint32_t crossfade_samples = 0) requires (
    PASS_TYPE == FilterPassType::LOW_PASS or PASS_TYPE == FilterPassType::HIGH_PASS);

    /**
     * @brief   Constructor (band-pass and band-stop).
     *
     * @param crossfade_samples  The number of samples to smooth the transition between old and new coefficients.
     */
    IIRCheby2(double normalized_lowcut_freq, double normalized_highcut_freq, double stopband_ripple_db,
              uint32_t crossfade_samples = 0) requires (
    PASS_TYPE == FilterPassType::BAND_PASS or PASS_TYPE == FilterPassType::BAND_STOP);

    /**
     * @brief   Configure the filter (low-pass and high-pass).
     *
     * @param normalized_cutoff_frequency  Normalized cutoff frequency.
     * @param stopband_ripple_db  Stopband ripple in dB.
     */
    void configure(double normalized_cutoff_frequency, double stopband_ripple_db) requires (
    PASS_TYPE == FilterPassType::LOW_PASS or PASS_TYPE == FilterPassType::HIGH_PASS);

    /**
     * @brief   Configure the filter (band-pass and band-stop).
     *
     * @param normalized_lowcut_freq  Normalized low-pass cutoff frequency.
     * @param normalized_highcut_freq  Normalized high-pass cutoff frequency.
     * @param stopband_ripple_db  Stopband ripple in dB.
     */
    void configure(double normalized_lowcut_freq, double normalized_highcut_freq, double stopband_ripple_db) requires (
    PASS_TYPE == FilterPassType::BAND_PASS or PASS_TYPE == FilterPassType::BAND_STOP);

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

    double _stopband_ripple_db = 0;
};


template<uint32_t N, typename T, FilterPassType PASS_TYPE>
IIRCheby2<N, T, PASS_TYPE>::IIRCheby2(double normalized_cutoff_frequency, double stopband_ripple_db,
                                      uint32_t crossfade_samples) requires (
PASS_TYPE == FilterPassType::LOW_PASS or PASS_TYPE == FilterPassType::HIGH_PASS)
        : IIRFilter<N, T, PASS_TYPE>(crossfade_samples) {
    configure(normalized_cutoff_frequency, stopband_ripple_db);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
IIRCheby2<N, T, PASS_TYPE>::IIRCheby2(double normalized_lowcut_freq, double normalized_highcut_freq,
                                      double stopband_ripple_db, uint32_t crossfade_samples) requires (
PASS_TYPE == FilterPassType::BAND_PASS or PASS_TYPE == FilterPassType::BAND_STOP)
        : IIRFilter<N, T, PASS_TYPE>(crossfade_samples) {
    configure(normalized_lowcut_freq, normalized_highcut_freq, stopband_ripple_db);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
void IIRCheby2<N, T, PASS_TYPE>::configure(double normalized_cutoff_frequency, double stopband_ripple_db) requires (
PASS_TYPE == FilterPassType::LOW_PASS or PASS_TYPE == FilterPassType::HIGH_PASS) {
    stopband_ripple_db = std::abs(stopband_ripple_db);
    if (stopband_ripple_db != _stopband_ripple_db) {
        _stopband_ripple_db = stopband_ripple_db;
        init_analog();
    }
    IIRFilter<N, T, PASS_TYPE>::calculate_cascades(normalized_cutoff_frequency);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
void IIRCheby2<N, T, PASS_TYPE>::configure(double normalized_lowcut_freq, double normalized_highcut_freq,
                                           double stopband_ripple_db) requires (
PASS_TYPE == FilterPassType::BAND_PASS or PASS_TYPE == FilterPassType::BAND_STOP) {
    stopband_ripple_db = std::abs(stopband_ripple_db);
    if (stopband_ripple_db != _stopband_ripple_db) {
        _stopband_ripple_db = stopband_ripple_db;
        init_analog();
    }
    IIRFilter<N, T, PASS_TYPE>::calculate_cascades(normalized_lowcut_freq, normalized_highcut_freq);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
double IIRCheby2<N, T, PASS_TYPE>::get_analog_gain() const {
    return 1.0;
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
void IIRCheby2<N, T, PASS_TYPE>::init_analog() {
    const double delta = 1.0 / std::sqrt(std::exp(_stopband_ripple_db * 0.1 * M_LN10) - 1);
    const double mu = std::asinh(1.0 / delta) / N;
    const double sinh_mu = std::sinh(mu);
    const double cosh_mu = std::cosh(mu);

    if constexpr (N & 1) {
        IIRFilter<N, T, PASS_TYPE>::_analog_pole_zero_pairs[(N + 1) / 2 - 1] = {
                -1.0 / sinh_mu,
                INFINITY_VALUE
        };
    }

    for (uint32_t i = 0; i < N / 2; ++i) {
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