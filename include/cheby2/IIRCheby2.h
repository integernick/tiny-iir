#pragma once

#include "common/IIRFilter.h"

namespace tiny_iir {
/**
 * @brief   IIR Chebyshev Type II filter.
 *
 * @tparam  N             The number of coefficients.
 * @tparam  T             The data type.
 * @tparam  PASS_TYPE     The filter pass type.
 * @tparam  DESIGN_T      The filter design type.
 */
template<uint32_t N = 2, typename T = double, FilterPassType PASS_TYPE = FilterPassType::LowPass,
        typename DESIGN_T = double>
class IIRCheby2 final : public IIRFilter<N, T, PASS_TYPE, DESIGN_T> {
    using DT = DESIGN_T;

public:
    /**
     * @brief   Constructor (low-pass and high-pass).
     *
     * @param normalized_cutoff_frequency   Normalized cutoff frequency.
     * @param stopband_ripple_db            Stopband attenuation in dB.
     * @param crossfade_samples             The number of samples to smooth the transition between old and new coefficients.
     */
    IIRCheby2(DT normalized_cutoff_frequency, DT stopband_ripple_db, uint32_t crossfade_samples = 0) requires (
    PASS_TYPE == FilterPassType::LowPass or PASS_TYPE == FilterPassType::HighPass);

    /**
     * @brief   Constructor (band-pass and band-stop).
     *
     * @param normalized_lowcut_freq    Normalized low-pass cutoff frequency.
     * @param normalized_highcut_freq   Normalized high-pass cutoff frequency.
     * @param stopband_ripple_db        Stopband attenuation in dB.
     * @param crossfade_samples         The number of samples to smooth the transition between old and new coefficients.
     */
    IIRCheby2(DT normalized_lowcut_freq, DT normalized_highcut_freq, DT stopband_ripple_db,
              uint32_t crossfade_samples = 0) requires (
    PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop);

    /**
     * @brief   Configure the filter (low-pass and high-pass).
     *
     * @param normalized_cutoff_frequency   Normalized cutoff frequency.
     * @param stopband_ripple_db            Stopband attenuation in dB.
     */
    void configure(DT normalized_cutoff_frequency, DT stopband_ripple_db) requires (
    PASS_TYPE == FilterPassType::LowPass or PASS_TYPE == FilterPassType::HighPass);

    /**
     * @brief   Configure the filter (band-pass and band-stop).
     *
     * @param normalized_lowcut_freq    Normalized low-pass cutoff frequency.
     * @param normalized_highcut_freq   Normalized high-pass cutoff frequency.
     * @param stopband_ripple_db        Stopband attenuation in dB.
     */
    void configure(DT normalized_lowcut_freq, DT normalized_highcut_freq, DT stopband_ripple_db) requires (
    PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop);

private:
    /**
     * @brief   Get the analog gain.
     *
     * @return  The analog gain.
     */
    [[nodiscard]] DT get_analog_gain() const final;

    /**
     * @brief   Initialize analog filter poles and zeros.
     */
    void init_analog() final;

    DT _stopband_ripple_db = 0;
};


template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
IIRCheby2<N, T, PASS_TYPE, DT>::IIRCheby2(DT normalized_cutoff_frequency, DT stopband_ripple_db,
                                          uint32_t crossfade_samples) requires (
PASS_TYPE == FilterPassType::LowPass or PASS_TYPE == FilterPassType::HighPass)
        : IIRFilter<N, T, PASS_TYPE, DT>(crossfade_samples) {
    configure(normalized_cutoff_frequency, stopband_ripple_db);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
IIRCheby2<N, T, PASS_TYPE, DT>::IIRCheby2(DT normalized_lowcut_freq, DT normalized_highcut_freq,
                                          DT stopband_ripple_db, uint32_t crossfade_samples) requires (
PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop)
        : IIRFilter<N, T, PASS_TYPE, DT>(crossfade_samples) {
    configure(normalized_lowcut_freq, normalized_highcut_freq, stopband_ripple_db);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRCheby2<N, T, PASS_TYPE, DT>::configure(DT normalized_cutoff_frequency, DT stopband_ripple_db) requires (
PASS_TYPE == FilterPassType::LowPass or PASS_TYPE == FilterPassType::HighPass) {
    stopband_ripple_db = std::abs(stopband_ripple_db);

    if (stopband_ripple_db != _stopband_ripple_db) {
        _stopband_ripple_db = stopband_ripple_db;
        init_analog();
    }

    IIRFilter<N, T, PASS_TYPE, DT>::calculate_cascades(normalized_cutoff_frequency);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRCheby2<N, T, PASS_TYPE, DT>::configure(DT normalized_lowcut_freq, DT normalized_highcut_freq,
                                               DT stopband_ripple_db) requires (
PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop) {
    stopband_ripple_db = std::abs(stopband_ripple_db);

    if (stopband_ripple_db != _stopband_ripple_db) {
        _stopband_ripple_db = stopband_ripple_db;
        init_analog();
    }

    IIRFilter<N, T, PASS_TYPE, DT>::calculate_cascades(normalized_lowcut_freq, normalized_highcut_freq);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
DT IIRCheby2<N, T, PASS_TYPE, DT>::get_analog_gain() const {
    return 1.0;
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRCheby2<N, T, PASS_TYPE, DT>::init_analog() {
    const DT delta = 1.0 / std::sqrt(std::exp(_stopband_ripple_db * 0.1 * std::numbers::ln10_v<DT>) - 1);
    const DT mu = std::asinh(1.0 / delta) / N;
    const DT sinh_mu = std::sinh(mu);
    const DT cosh_mu = std::cosh(mu);

    if constexpr (N & 1) {
        IIRFilter<N, T, PASS_TYPE, DT>::_analog_pole_zero_pairs[(N + 1) / 2 - 1] = {
                -DT{1} / sinh_mu,
                std::numeric_limits<DT>::infinity()
        };
    }

    for (uint32_t i = 0; i < N / 2; ++i) {
        // Angle from the imaginary axis:
        const DT phi = static_cast<DT>(2 * i + 1) * std::numbers::pi_v<DT> * static_cast<DT>(0.5) / N;
        const DT sin_phi = std::sin(phi);
        const DT cos_phi = std::cos(phi);

        const DT pole_s_real = -sinh_mu * sin_phi;
        const DT pole_s_imag = cosh_mu * cos_phi;
        const DT pole_squared = pole_s_real * pole_s_real + pole_s_imag * pole_s_imag;
        const DT zero_imag = DT{1} / cos_phi;

        IIRFilter<N, T, PASS_TYPE, DT>::_analog_pole_zero_pairs[i] = {
                Complex{pole_s_real / pole_squared, pole_s_imag / pole_squared},
                Complex{DT{0}, zero_imag}
        };
    }
}

}