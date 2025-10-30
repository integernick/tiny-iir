#pragma once

#include "common/IIRFilter.h"

namespace tiny_iir {

/**
 * @brief   IIR Chebyshev Type I filter.
 *
 * @tparam  N             The number of coefficients.
 * @tparam  T             The data type.
 * @tparam  PASS_TYPE     The filter pass type.
 * @tparam  DESIGN_T      The filter design type.
 */
template<uint32_t N = 2, typename T = double, FilterPassType PASS_TYPE = FilterPassType::LowPass,
        typename DESIGN_T = double>
class IIRCheby1 final : public IIRFilter<N, T, PASS_TYPE, DESIGN_T> {
    using DT = DESIGN_T;

public:
    /**
     * @brief   Constructor (low-pass and high-pass).
     *
     * @param normalized_cutoff_frequency   Normalized cutoff frequency.
     * @param passband_ripple_db            Passband ripple in dB.
     * @param crossfade_samples             The number of samples to smooth the transition between old and new coefficients.
     */
    IIRCheby1(DT normalized_cutoff_frequency, DT passband_ripple_db, uint32_t crossfade_samples = 0) requires (
    PASS_TYPE == FilterPassType::LowPass or PASS_TYPE == FilterPassType::HighPass);

    /**
     * @brief   Constructor (band-pass and band-stop).
     *
     * @param normalized_lowcut_freq    Normalized low-pass cutoff frequency.
     * @param normalized_highcut_freq   Normalized high-pass cutoff frequency.
     * @param passband_ripple_db        Passband ripple in dB.
     * @param crossfade_samples         The number of samples to smooth the transition between old and new coefficients.
     */
    IIRCheby1(DT normalized_lowcut_freq, DT normalized_highcut_freq, DT passband_ripple_db,
              uint32_t crossfade_samples = 0) requires (
    PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop);

    /**
     * @brief   Configure the filter (low-pass and high-pass).
     *
     * @param normalized_cutoff_frequency   Normalized cutoff frequency.
     * @param passband_ripple_db            Passband ripple in dB.
     */
    void configure(DT normalized_cutoff_frequency, DT passband_ripple_db) requires (
    PASS_TYPE == FilterPassType::LowPass or PASS_TYPE == FilterPassType::HighPass);

    /**
     * @brief   Configure the filter (band-pass and band-stop).
     *
     * @param normalized_lowcut_freq    Normalized low-pass cutoff frequency.
     * @param normalized_highcut_freq   Normalized high-pass cutoff frequency.
     * @param passband_ripple_db        Passband ripple in dB.
     */
    void configure(DT normalized_lowcut_freq, DT normalized_highcut_freq, DT passband_ripple_db) requires (
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

    DT _passband_ripple_db = 0;
};


template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
IIRCheby1<N, T, PASS_TYPE, DT>::IIRCheby1(DT normalized_cutoff_frequency, DT passband_ripple_db,
                                          uint32_t crossfade_samples) requires (
PASS_TYPE == FilterPassType::LowPass or PASS_TYPE == FilterPassType::HighPass)
        : IIRFilter<N, T, PASS_TYPE, DT>(crossfade_samples) {
    configure(normalized_cutoff_frequency, passband_ripple_db);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
IIRCheby1<N, T, PASS_TYPE, DT>::IIRCheby1(DT normalized_lowcut_freq, DT normalized_highcut_freq,
                                          DT passband_ripple_db, uint32_t crossfade_samples) requires (
PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop)
        : IIRFilter<N, T, PASS_TYPE, DT>(crossfade_samples) {
    configure(normalized_lowcut_freq, normalized_highcut_freq, passband_ripple_db);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRCheby1<N, T, PASS_TYPE, DT>::configure(DT normalized_cutoff_frequency, DT passband_ripple_db) requires (
PASS_TYPE == FilterPassType::LowPass or PASS_TYPE == FilterPassType::HighPass) {
    passband_ripple_db = std::abs(passband_ripple_db);

    if (passband_ripple_db != _passband_ripple_db) {
        _passband_ripple_db = passband_ripple_db;
        init_analog();
    }

    IIRFilter<N, T, PASS_TYPE, DT>::calculate_cascades(normalized_cutoff_frequency);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRCheby1<N, T, PASS_TYPE, DT>::configure(DT normalized_lowcut_freq, DT normalized_highcut_freq,
                                               DT passband_ripple_db) requires (
PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop) {
    passband_ripple_db = std::abs(passband_ripple_db);

    if (passband_ripple_db != _passband_ripple_db) {
        _passband_ripple_db = passband_ripple_db;
        init_analog();
    }

    IIRFilter<N, T, PASS_TYPE, DT>::calculate_cascades(normalized_lowcut_freq, normalized_highcut_freq);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
DT IIRCheby1<N, T, PASS_TYPE, DT>::get_analog_gain() const {
    if constexpr (N & 1) {
        return 1.0;
    } else {
        return std::exp(-_passband_ripple_db / 20 * M_LN10);
    }
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRCheby1<N, T, PASS_TYPE, DT>::init_analog() {
    const DT epsilon = std::sqrt(std::exp(_passband_ripple_db * 0.1 * M_LN10) - 1);
    const DT mu = std::asinh(DT{1} / epsilon) / N;
    const DT sinh_mu = std::sinh(mu);
    const DT cosh_mu = std::cosh(mu);

    if constexpr (N & 1) {
        IIRFilter<N, T, PASS_TYPE, DT>::_analog_pole_zero_pairs[(N + 1) / 2 - 1] = {
                Complex<DT>{-sinh_mu, DT{0}},
                Complex<DT>{std::numeric_limits<DT>::infinity(), DT{0}}
        };
    }

    for (uint32_t i = 0; i < N / 2; ++i) {
        // Angle from the imaginary axis
        const DT phi = static_cast<DT>(2 * i + 1) * std::numbers::pi_v<DT> * static_cast<DT>(0.5) / N;
        const DT sin_phi = std::sin(phi);
        const DT cos_phi = std::cos(phi);

        const DT pole_s_real = -sinh_mu * sin_phi;
        const DT pole_s_imag = cosh_mu * cos_phi;
        IIRFilter<N, T, PASS_TYPE, DT>::_analog_pole_zero_pairs[i] = {
                Complex{pole_s_real, pole_s_imag},
                Complex{std::numeric_limits<DT>::infinity(), DT{0}}
        };
    }
}

}