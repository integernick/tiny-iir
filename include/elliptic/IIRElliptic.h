#pragma once

#include "elliptic_utils.h"
#include "common/IIRFilter.h"

namespace tiny_iir {

/**
 * @brief Elliptic filter
 *
 * @tparam N   Filter order.
 * @tparam T   Data type.
 * @tparam PASS_TYPE   Pass type (low, high).
 * @tparam DESIGN_T    The filter design type.
 */
template<uint32_t N = 2, typename T = double, FilterPassType PASS_TYPE = FilterPassType::LowPass,
        typename DESIGN_T = double>
class IIRElliptic : public IIRFilter<N, T, PASS_TYPE, DESIGN_T> {
    static_assert(std::is_same_v<DESIGN_T, float> or std::is_same_v<DESIGN_T, double>,
                  "DESIGN_T must be float or double");
    using DT = DESIGN_T;

    using Complex = Complex<DT>;

public:
    /**
     * @brief   Constructor (low-pass and high-pass).
     *
     * @param normalized_cutoff_frequency   Normalized cutoff frequency.
     * @param pass_ripple_db                Passband ripple in dB.
     * @param stop_attenuation_db           Stopband attenuation in dB.
     * @param crossfade_samples             The number of samples to smooth the transition between old and new coefficients.
     */
    IIRElliptic(DT normalized_cutoff_frequency, DT pass_ripple_db, DT stop_attenuation_db,
                uint32_t crossfade_samples = 0) requires (
    PASS_TYPE == FilterPassType::LowPass or PASS_TYPE == FilterPassType::HighPass);

    /**
     * @brief   Constructor (band-pass and band-stop).
     *
     * @param normalized_lowcut_freq    Normalized low-pass cutoff frequency.
     * @param normalized_highcut_freq   Normalized high-pass cutoff frequency.
     * @param pass_ripple_db            Passband ripple in dB.
     * @param stop_attenuation_db       Stopband attenuation in dB.
     * @param crossfade_samples         The number of samples to smooth the transition between old and new coefficients.
     */
    IIRElliptic(DT normalized_lowcut_freq, DT normalized_highcut_freq,
                DT pass_ripple_db, DT stop_attenuation_db, uint32_t crossfade_samples = 0) requires (
    PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop);

    /**
     * @brief   Configure the filter (low-pass and high-pass).
     *
     * @param normalized_cutoff_frequency   Normalized cutoff frequency.
     * @param pass_ripple_db                Passband ripple in dB.
     * @param stop_attenuation_db           Stopband attenuation in dB.
     */
    void configure(DT normalized_cutoff_frequency, DT pass_ripple_db, DT stop_attenuation_db) requires (
    PASS_TYPE == FilterPassType::LowPass or PASS_TYPE == FilterPassType::HighPass);

    /**
     * @brief   Configure the filter (band-pass and band-stop).
     *
     * @param normalized_lowcut_freq    Normalized low-pass cutoff frequency.
     * @param normalized_highcut_freq   Normalized high-pass cutoff frequency.
     * @param pass_ripple_db            Passband ripple in dB.
     * @param stop_attenuation_db       Stopband attenuation in dB.
     */
    void configure(DT normalized_lowcut_freq, DT normalized_highcut_freq,
                   DT pass_ripple_db, DT stop_attenuation_db) requires (
    PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop);

private:
    static constexpr uint32_t L = N / 2;
    static constexpr Complex IMAG_UNIT = {DT{0}, DT{1}};

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

    DT _pass_ripple_db = 0;
    DT _stop_attenuation_db = 0;
};


template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
IIRElliptic<N, T, PASS_TYPE, DT>::IIRElliptic(DT normalized_cutoff_frequency, DT pass_ripple_db,
                                              DT stop_attenuation_db, uint32_t crossfade_samples) requires (
PASS_TYPE == FilterPassType::LowPass or PASS_TYPE == FilterPassType::HighPass)
        : IIRFilter<N, T, PASS_TYPE, DT>(crossfade_samples) {
    configure(normalized_cutoff_frequency, pass_ripple_db, stop_attenuation_db);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
IIRElliptic<N, T, PASS_TYPE, DT>::IIRElliptic(DT normalized_lowcut_freq, DT normalized_highcut_freq,
                                              DT pass_ripple_db, DT stop_attenuation_db,
                                              uint32_t crossfade_samples) requires (
PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop)
        : IIRFilter<N, T, PASS_TYPE, DT>(crossfade_samples) {
    configure(normalized_lowcut_freq, normalized_highcut_freq, pass_ripple_db, stop_attenuation_db);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRElliptic<N, T, PASS_TYPE, DT>::configure(DT normalized_cutoff_frequency,
                                                 DT pass_ripple_db, DT stop_attenuation_db) requires (
PASS_TYPE == FilterPassType::LowPass or PASS_TYPE == FilterPassType::HighPass) {
    pass_ripple_db = std::abs(pass_ripple_db);
    stop_attenuation_db = std::abs(stop_attenuation_db);

    if (pass_ripple_db != _pass_ripple_db || stop_attenuation_db != _stop_attenuation_db) {
        _pass_ripple_db = pass_ripple_db;
        _stop_attenuation_db = stop_attenuation_db;
        init_analog();
    }

    IIRFilter<N, T, PASS_TYPE, DT>::calculate_cascades(normalized_cutoff_frequency);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRElliptic<N, T, PASS_TYPE, DT>::configure(DT normalized_lowcut_freq, DT normalized_highcut_freq,
                                                 DT pass_ripple_db, DT stop_attenuation_db) requires (
PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop) {
    pass_ripple_db = std::abs(pass_ripple_db);
    stop_attenuation_db = std::abs(stop_attenuation_db);

    if (pass_ripple_db != _pass_ripple_db || stop_attenuation_db != _stop_attenuation_db) {
        _pass_ripple_db = pass_ripple_db;
        _stop_attenuation_db = stop_attenuation_db;
        init_analog();
    }

    IIRFilter<N, T, PASS_TYPE, DT>::calculate_cascades(normalized_lowcut_freq, normalized_highcut_freq);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
DT IIRElliptic<N, T, PASS_TYPE, DT>::get_analog_gain() const {
    if constexpr (N & 1) {
        return DT{1};
    } else {
        return std::exp(-_pass_ripple_db / DT{20} * std::numbers::ln10_v<DT>);
    }
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRElliptic<N, T, PASS_TYPE, DT>::init_analog() {
    const DT eps_p = std::sqrt(std::exp(_pass_ripple_db * DT{0.1} * std::numbers::ln10_v<DT>) - DT{1});
    const DT eps_s = std::sqrt(std::exp(_stop_attenuation_db * DT{0.1} * std::numbers::ln10_v<DT>) - DT{1});

    const DT k1 = eps_p / eps_s; // The ratio (eps_p / eps_s) << 1

    // Solution of sn(j*v0*N*K1, k1) = j/eps_p
    const DT K1 = calculate_elliptic_integral(k1);
    const DT K1_prime = calculate_elliptic_integral(get_complimentary(k1));
    const DT R1 = K1_prime / K1;
    const Complex v0_c = -IMAG_UNIT * asn<DT>(IMAG_UNIT / eps_p, k1, R1) / static_cast<DT>(N);
    const DT v0 = v0_c.real();

    const DT k = solve_degree_equation<DT>(N, get_complimentary(k1)); // The ratio (w_p / w_s) < 1

    if (k > DT{1}) {
        return;
    }

    if constexpr (N & 1) {
        const Complex pole = IMAG_UNIT * sn<DT>(IMAG_UNIT * v0, k);
        const Complex zero = Complex{std::numeric_limits<DT>::infinity(), DT{0}};
        IIRFilter<N, T, PASS_TYPE, DT>::_analog_pole_zero_pairs[(N + 1) / 2 - 1] = {pole, zero};
    }

    for (uint32_t i = 0; i < N / 2; ++i) {
        const DT u_i = (DT{2} * i + DT{1}) / N;
        const DT zeta = cd<DT>(u_i, k).real();
        const Complex zero = Complex{DT{0}, DT{1} / (k * zeta)};
        const Complex pole = IMAG_UNIT * cd<DT>(u_i - IMAG_UNIT * v0, k);

        IIRElliptic<N, T, PASS_TYPE, DT>::_analog_pole_zero_pairs[i] = {pole, zero};
    }
}

}