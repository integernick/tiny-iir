#pragma once

#include "CascadeFilter.h"
#include "PassTypeData.h"

namespace tiny_iir {

/**
 * @brief   IIR filter.
 *
 * @tparam N   Analog prototype filter order.
 * @tparam T   Data type.
 * @tparam PASS_TYPE   Pass type (low, high, band-pass, band-stop).
 * @tparam DESIGN_T    The filter design type.
 */
template<uint32_t N, typename T = double, FilterPassType PASS_TYPE = FilterPassType::LowPass,
        typename DESIGN_T = double>
class IIRFilter {
    static_assert(std::is_same_v<T, double> or std::is_same_v<T, float>
                  or std::is_same_v<T, q31_t> or std::is_same_v<T, q15_t>,
                  "T must be double, float, q31_t, or q15_t");

    static_assert(std::is_same_v<DESIGN_T, double> or std::is_same_v<DESIGN_T, float>,
                  "DESIGN_T must be float or double");

    using DT = DESIGN_T;

public:
    using ValueType = T;

    /**
     * @brief   The destructor.
     */
    virtual ~IIRFilter() = default;

    /**
     * @brief   Process a single sample.
     *
     * @tparam U  The type of the input sample.
     * @param x The input sample.
     * @return  The output sample.
     */
    template<typename U>
    [[nodiscard]] U process(U x);

    /**
     * @brief   Process a batch of samples.
     *
     * @param[in] x             The input buffer.
     * @param[out] out          The output buffer.
     * @paramp[in] num_samples  The number of samples to process.
     */
    void process(const T *x, T *out, uint32_t num_samples);

    /**
     * @brief   Process a batch of samples.
     *
     * @param[in] x         The input buffer.
     * @param num_samples   The number of samples to process.
     * @return  The final output sample.
     */
    [[nodiscard]] T process(const T *x, uint32_t num_samples);

    /**
     * @brief   Process a batch of samples (input is DT, native is not DT).
     *
     * @param[in] x             The input buffer.
     * @param[out] out          The output buffer.
     * @paramp[in] num_samples  The number of samples to process.
     */
    template<typename U>
    void process(const U *x, T *out, uint32_t num_samples) requires (
    std::is_same_v<U, DT> and !std::is_same_v<T, DT>);

    /**
     * @brief   Process a batch of samples (input is DT, native is not DT).
     *
     * @param[in] x         The input buffer.
     * @param num_samples   The number of samples to process.
     * @return  The final output sample.
     */
    template<typename U>
    [[nodiscard]] U process(const U *x, uint32_t num_samples) requires (
    std::is_same_v<U, DT> and !std::is_same_v<T, DT>);

    /**
     * @brief   Get the filter gain.
     *
     * @return  The filter gain.
     */
    [[nodiscard]] DT get_gain() const;

    /**
     * @brief   Get filter coefficients.
     *
     * @return  Pointer to the coefficients array.
     */
    [[nodiscard]] const T *get_coefficients() const;

    /**
     * @brief   Get the biquad coefficients (DT representation).
     *
     * @param biquad_idx  The index of the biquad.
     * @return  The biquad coefficients.
     */
    [[nodiscard]] BiquadCoefficients<DT> get_biquad_coefficients(uint32_t biquad_idx) const;

    /**
     * @brief   Print filter coefficients.
     */
    void print_coefficients() const;

    /**
     * @brief   Set the number of samples to smooth the transition between old and new coefficients.
     *
     * @param crossfade_samples  The number of samples to smooth the transition between old and new coefficients.
     */
    void set_crossfade_samples(uint32_t crossfade_samples);

    /**
     * @brief   Reset filter state.
     */
    void reset_state();

    static constexpr uint32_t REAL_ORDER = PassTypeData<PASS_TYPE, N, DT>::CASCADE_ORDER;

    // Not just CascadeFilter<N, T>::NUMBER_OF_BIQUAD_BLOCKS, because the order is doubled for band-pass/band-stop
    static constexpr uint32_t NUMBER_OF_BIQUAD_BLOCKS
            = CascadeFilter<PassTypeData<PASS_TYPE, N, DT>::CASCADE_ORDER, T, DT>::NUMBER_OF_BIQUAD_BLOCKS;

protected:
    /**
     * @brief   Initialize analog filter poles and zeros.
     */
    virtual void init_analog() = 0;

    /**
     * @brief   Get the analog gain.
     *
     * @return  The analog gain.
     */
    [[nodiscard]] virtual DT get_analog_gain() const = 0;

    /**
     * @brief   Constructor.
     */
    explicit IIRFilter(uint32_t crossfade_samples = 0);

    /**
     * @brief   Calculate biquad cascades coefficients (low-pass and high-pass).
     *
     * @param normalized_cutoff_frequency  Normalized cutoff frequency.
     */
    void calculate_cascades(DT normalized_cutoff_frequency) requires (
    PASS_TYPE == FilterPassType::LowPass or PASS_TYPE == FilterPassType::HighPass);

    /**
     * @brief   Calculate biquad cascades coefficients (band-pass and band-stop).
     *
     * @param cutlow_freq  Low cutoff frequency.
     * @param cuthigh_freq  High cutoff frequency.
     */
    void calculate_cascades(DT cutlow_freq, DT cuthigh_freq) requires (
    PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop);

    CascadeFilter<PassTypeData<PASS_TYPE, N, DT>::CASCADE_ORDER, T, DT> _cascade_filter;

    // Pole/zero pairs without the conjugates (last one is the pair with pole on the real axis)
    PoleZeroPair<DT> _analog_pole_zero_pairs[(N + 1) / 2];

private:
    struct FrequencyConfig {
        DT w1{};
        DT w2{};
    };

    /**
     * @brief   Add a pair of S-plane pole-zero conjugates.
     *
     * @param pole_zero_pair  The pole-zero conjugates pair.
     */
    void add_pole_zero_conjugates_pair(const PoleZeroPair<DT> &pole_zero_pair);

    /**
     * @brief   Add a single S-plane pole-zero pair (pole and zero lie on the real axis).
     *
     * @param pole_zero_pair  The pole-zero pair.
     */
    void add_pole_zero_single_pair(const PoleZeroPair<DT> &pole_zero_pair);

    /**
     * @brief   Push a new block of biquad coefficients.
     *
     * @param biquad_coefficients  The biquad coefficients.
     * @return  True if pushing the coefficients was successful, false otherwise.
     */
    void push_biquad_coefficients(BiquadCoefficients<DT> &biquad_coefficients);

    /**
     * @brief   Add a pair of S-plane pole-zero pairs (band-pass/band-stop filters with an odd order).
     *
     * @note    This one is needed for band-pass/band-stop filters with an odd design order;
     *          Each real-axis pole (zero) transforms to two poles (zeros) in the complex plane after the transform,
     *          after which the poles and zeros are arranged in pairs.
     * @param pole_zero_pairs  The pole-zero pairs.
     */
    void add_pole_zero_pairs(const std::pair<PoleZeroPair<DT>, PoleZeroPair<DT>> &pole_zero_pairs) requires (
    (PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop) and (N & 1) != 0);

    PassTypeData<PASS_TYPE, N, DT> _pass_type_data;
    FrequencyConfig _frequency_config;
};

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
DT IIRFilter<N, T, PASS_TYPE, DT>::get_gain() const {
    return _cascade_filter.get_gain();
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
const T *IIRFilter<N, T, PASS_TYPE, DT>::get_coefficients() const {
    return _cascade_filter.get_coefficients();
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
BiquadCoefficients<DT> IIRFilter<N, T, PASS_TYPE, DT>::get_biquad_coefficients(uint32_t biquad_idx) const {
    return _cascade_filter.get_biquad_coefficients(biquad_idx);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRFilter<N, T, PASS_TYPE, DT>::print_coefficients() const {
    _cascade_filter.print_coefficients();
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
T IIRFilter<N, T, PASS_TYPE, DT>::process(const T *x, uint32_t num_samples) {
    return _cascade_filter.process(x, num_samples);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
template<typename U>
void IIRFilter<N, T, PASS_TYPE, DT>::process(const U *x, T *out, uint32_t num_samples) requires (
std::is_same_v<U, DT> and !std::is_same_v<T, DT>) {
    _cascade_filter.process(x, out, num_samples);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
template<typename U>
U IIRFilter<N, T, PASS_TYPE, DT>::process(const U *x, uint32_t num_samples) requires (
std::is_same_v<U, DT> and !std::is_same_v<T, DT>) {
    return _cascade_filter.process(x, num_samples);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRFilter<N, T, PASS_TYPE, DT>::set_crossfade_samples(uint32_t crossfade_samples) {
    _cascade_filter.set_crossfade_samples(crossfade_samples);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRFilter<N, T, PASS_TYPE, DT>::reset_state() {
    _cascade_filter.reset_state();
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
IIRFilter<N, T, PASS_TYPE, DT>::IIRFilter(uint32_t crossfade_samples) {
    set_crossfade_samples(crossfade_samples);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRFilter<N, T, PASS_TYPE, DT>::calculate_cascades(DT normalized_cutoff_frequency) requires (
PASS_TYPE == FilterPassType::LowPass or PASS_TYPE == FilterPassType::HighPass) {
    if (normalized_cutoff_frequency == _frequency_config.w1) {
        return;
    }

    _frequency_config.w1 = normalized_cutoff_frequency;

    _pass_type_data.init(normalized_cutoff_frequency);

    _cascade_filter.set_gain(get_analog_gain());
    _cascade_filter.reset();
    _cascade_filter.init_coefficients();

    // Start from the lowest Q factor pole (closest to the real axis)
    if constexpr (N & 1) {
        PoleZeroPair pole_zero_pair = _analog_pole_zero_pairs[(N + 1) / 2 - 1];
        pole_zero_pair.pole = _pass_type_data.transform(pole_zero_pair.pole);
        pole_zero_pair.zero = _pass_type_data.transform(pole_zero_pair.zero);
        add_pole_zero_single_pair(pole_zero_pair);
    }

    // Start from the lowest Q factor pole (closest to the real axis)
    for (int i = static_cast<int>(N / 2) - 1; i >= 0; --i) {
        PoleZeroPair pole_zero_pair = _analog_pole_zero_pairs[i];
        pole_zero_pair.pole = _pass_type_data.transform(pole_zero_pair.pole);
        pole_zero_pair.zero = _pass_type_data.transform(pole_zero_pair.zero);
        add_pole_zero_conjugates_pair(pole_zero_pair);
    }

    _cascade_filter.init();
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRFilter<N, T, PASS_TYPE, DT>::calculate_cascades(DT cutlow_freq, DT cuthigh_freq) requires (
PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop) {
    if (cuthigh_freq == _frequency_config.w2 && cutlow_freq == _frequency_config.w1) {
        return;
    }

    _frequency_config.w1 = cutlow_freq;
    _frequency_config.w2 = cuthigh_freq;

    _cascade_filter.set_gain(get_analog_gain());
    _cascade_filter.reset_blocks();
    _cascade_filter.reset();
    _cascade_filter.init_coefficients();

    constexpr DT MIN_FREQ = DT{1e-8};
    constexpr DT MAX_FREQ = DT{2} * std::numbers::pi_v<DT> - MIN_FREQ;

    if (cutlow_freq < MIN_FREQ) {
        cutlow_freq = MIN_FREQ;
    }

    if (cuthigh_freq > MAX_FREQ) {
        cuthigh_freq = MAX_FREQ;
    }

    _pass_type_data.init(cutlow_freq, cuthigh_freq);

    auto add_two_pole_zero_pairs_from_one = [&](const PoleZeroPair<DT> &pole_zero_pair) {
        auto [pole1, pole2] = _pass_type_data.transform(pole_zero_pair.pole);
        auto [zero1, zero2] = _pass_type_data.transform(pole_zero_pair.zero);

        PoleZeroPair pz1(pole1, zero1);
        add_pole_zero_conjugates_pair(pz1);

        PoleZeroPair pz2(pole2, zero2);
        add_pole_zero_conjugates_pair(pz2);
    };

    if constexpr (N & 1) {
        const PoleZeroPair<DT> &pole_zero_pair = _analog_pole_zero_pairs[(N + 1) / 2 - 1];
        auto [pole1, pole2] = _pass_type_data.transform(pole_zero_pair.pole);
        auto [zero1, zero2] = _pass_type_data.transform(pole_zero_pair.zero);

        PoleZeroPair pz1(pole1, zero1);
        PoleZeroPair pz2(pole2, zero2);
        add_pole_zero_pairs({pz1, pz2});
    }

    // Start from the lowest Q factor pole (closest to the real axis)
    for (int i = static_cast<int>(N / 2) - 1; i >= 0; --i) {
        add_two_pole_zero_pairs_from_one(_analog_pole_zero_pairs[i]);
    }

    _cascade_filter.init();
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRFilter<N, T, PASS_TYPE, DT>::add_pole_zero_conjugates_pair(const PoleZeroPair<DT> &pole_zero_pair) {
    const Complex<DT> &pole_z = pole_zero_pair.pole;
    const Complex<DT> &zero_z = pole_zero_pair.zero;

    // (z - z0)(z - z0*)=z^2 - z(z0+z0*) + |z0|^2 = z^2 * (1 - 2*Re(z0) * z^-1 + |z0|^2 * z^-2)
    BiquadCoefficients biquad_coefficients = {
            .b0 = DT{1},
            .b1 = -DT{2} * zero_z.real(),
            .b2 = zero_z.real() * zero_z.real() + zero_z.imag() * zero_z.imag(),
            .a1 = -DT{2} * pole_z.real(),
            .a2 = pole_z.real() * pole_z.real() + pole_z.imag() * pole_z.imag()
    };

    push_biquad_coefficients(biquad_coefficients);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRFilter<N, T, PASS_TYPE, DT>::add_pole_zero_single_pair(const PoleZeroPair<DT> &pole_zero_pair) {
    const Complex<DT> &pole_z = pole_zero_pair.pole;
    const Complex<DT> &zero_z = pole_zero_pair.zero;

    // (z - z0) = z * (1 - z^-1 * z0 + 0)
    BiquadCoefficients biquad_coefficients = {
            .b0 = DT{1},
            .b1 = -zero_z.real(),
            .b2 = DT{0},
            .a1 = -pole_z.real(),
            .a2 = DT{0}
    };

    push_biquad_coefficients(biquad_coefficients);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void IIRFilter<N, T, PASS_TYPE, DT>::push_biquad_coefficients(BiquadCoefficients<DT> &biquad_coefficients) {
    _cascade_filter.push_biquad_coefficients(biquad_coefficients);
    const DT biquad_gain = _pass_type_data.calculate_gain(biquad_coefficients);
    _cascade_filter.update_biquad_gain(DT{1} / biquad_gain);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
void
IIRFilter<N, T, PASS_TYPE, DT>::add_pole_zero_pairs(
        const std::pair<PoleZeroPair<DT>, PoleZeroPair<DT>> &pole_zero_pairs) requires (
(PASS_TYPE == FilterPassType::BandPass or PASS_TYPE == FilterPassType::BandStop) and (N & 1) != 0) {
    // (z - z1)(z - z2)=z^2 - z(z1+z2) + z1z2 = z^2 * (1 - (z1+z2)z^-1 + (z1*z2)*z^-2)
    const Complex<DT> &p1 = pole_zero_pairs.first.pole;
    const Complex<DT> &p2 = pole_zero_pairs.second.pole;
    const Complex<DT> &z1 = pole_zero_pairs.first.zero;
    const Complex<DT> &z2 = pole_zero_pairs.second.zero;

    BiquadCoefficients biquad_coefficients = {
            .b0 = DT{1},
            .b1 = -(pole_zero_pairs.first.zero.real() + pole_zero_pairs.second.zero.real()),
            .b2 = z1.imag() == DT{0}
                  ? z1.real() * z2.real()
                  : std::abs(z1),
            .a1 = -(pole_zero_pairs.first.pole.real() + pole_zero_pairs.second.pole.real()),
            .a2 = p1.imag() == DT{0}
                  ? p1.real() * p2.real()
                  : std::norm(p1)
    };

    push_biquad_coefficients(biquad_coefficients);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE, typename DT>
template<typename U>
U IIRFilter<N, T, PASS_TYPE, DT>::process(U x) {
    _cascade_filter.update_coefficients_crossfade();
    return _cascade_filter.process(x);
}
}