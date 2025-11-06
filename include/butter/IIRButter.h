#pragma once

#include "common/IIRFilter.h"

namespace tiny_iir {

template<uint32_t N = 2, typename T = double, FilterPassType PASS_TYPE = FilterPassType::LowPass,
        typename DESIGN_T = double>
class IIRButter final : public IIRFilter<N, T, PASS_TYPE, DESIGN_T> {
    using DT = DESIGN_T;

public:
    using IIRFilterBase = IIRFilter<N, T, PASS_TYPE, DT>;

    /**
     * @brief   Constructor (low-pass and high-pass).
     *
     * @param normalized_cutoff_frequency   Normalized cutoff frequency.
     * @param crossfade_samples             The number of samples to smooth the transition between old and new coefficients.
     */
    explicit IIRButter(DT normalized_cutoff_frequency, uint32_t crossfade_samples = 0) requires (
    PASS_TYPE == FilterPassType::LowPass or PASS_TYPE == FilterPassType::HighPass);

    /**
     * @brief   Constructor (band-pass and band-stop).
     *
     * @param normalized_lowcut_freq    Normalized low-pass cutoff frequency.
     * @param normalized_highcut_freq   Normalized high-pass cutoff frequency.
     * @param crossfade_samples         The number of samples to smooth the transition between old and new coefficients.
     */
    IIRButter(DT normalized_lowcut_freq, DT normalized_highcut_freq, uint32_t crossfade_samples = 0) requires (
    PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop);

    /**
     * @brief   Configure the filter (low-pass and high-pass).
     *
     * @param normalized_cutoff_frequency  Normalized cutoff frequency.
     */
    void configure(DT normalized_cutoff_frequency) requires (
    PASS_TYPE == FilterPassType::LowPass or PASS_TYPE == FilterPassType::HighPass) {
        init_analog();
        IIRFilterBase::calculate_cascades(normalized_cutoff_frequency);
    }

    /**
     * @brief   Configure the filter (band-pass and band-stop).
     *
     * @param normalized_lowcut_freq    Normalized low-pass cutoff frequency.
     * @param normalized_highcut_freq   Normalized high-pass cutoff frequency.
     */
    void configure(DT normalized_lowcut_freq, DT normalized_highcut_freq) requires (
    PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop) {
        init_analog();
        IIRFilterBase::calculate_cascades(normalized_lowcut_freq, normalized_highcut_freq);
    }

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
};


template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
IIRButter<N, T, PASS_TYPE, DT>::IIRButter(DT normalized_cutoff_frequency, uint32_t crossfade_samples) requires (
PASS_TYPE == FilterPassType::LowPass or PASS_TYPE == FilterPassType::HighPass)
        : IIRFilterBase(crossfade_samples) {
    configure(normalized_cutoff_frequency);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
IIRButter<N, T, PASS_TYPE, DT>::IIRButter(DT normalized_lowcut_freq, DT normalized_highcut_freq,
                                          uint32_t crossfade_samples) requires (
PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop)
        : IIRFilterBase(crossfade_samples) {
    configure(normalized_lowcut_freq, normalized_highcut_freq);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
DT IIRButter<N, T, PASS_TYPE, DT>::get_analog_gain() const {
    return DT{1};
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRButter<N, T, PASS_TYPE, DT>::init_analog() {
    if constexpr (N & 1) {
        IIRFilterBase::_analog_pole_zero_pairs[(N + 1) / 2 - 1] = {
                -DT{1},
                std::numeric_limits<DT>::infinity()
        };
    }

    for (uint32_t i = 0; i < N / 2; ++i) {
        // Angle from the imaginary axis:
        const DT phi = static_cast<DT>(2 * i + 1) * (std::numbers::pi_v<DT> / DT{2}) / N;
        const DT pole_s_real = -std::sin(phi);
        const DT pole_s_imag = std::cos(phi);
        IIRFilterBase::_analog_pole_zero_pairs[i] = {
                {pole_s_real,                         pole_s_imag},
                {std::numeric_limits<DT>::infinity(), DT{0}}
        };
    }
}

}