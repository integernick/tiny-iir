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
 */
template<uint32_t N,
        typename T = double,
        FilterPassType PASS_TYPE = FilterPassType::LOW_PASS>
class IIRFilter {
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
    [[nodiscard]] T process(U x);

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
     * @brief   Process a batch of samples (input is double, native is not double).
     *
     * @param[in] x             The input buffer.
     * @param[out] out          The output buffer.
     * @paramp[in] num_samples  The number of samples to process.
     */
    template<typename U, typename = std::enable_if_t<(std::is_same_v<U, double>
                                                      && !std::is_same_v<T, double>)>>
    void process(const U *x, T *out, uint32_t num_samples);

    /**
     * @brief   Process a batch of samples (input is double, native is not double).
     *
     * @param[in] x         The input buffer.
     * @param num_samples   The number of samples to process.
     * @return  The final output sample.
     */
    template<typename U, typename = std::enable_if_t<(std::is_same_v<U, double>
                                                      && !std::is_same_v<T, double>)>>
    [[nodiscard]] T process(const U *x, uint32_t num_samples);

    /**
     * @brief   Get the filter gain.
     *
     * @return  The filter gain.
     */
    [[nodiscard]] T get_gain() const;

    /**
     * @brief   Get filter coefficients.
     *
     * @return  Pointer to the coefficients array.
     */
    [[nodiscard]] const T *get_coefficients() const;

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

    static constexpr uint32_t REAL_ORDER = PassTypeData<PASS_TYPE, N>::CASCADE_ORDER;

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
    [[nodiscard]] virtual double get_analog_gain() const = 0;

    /**
     * @brief   Constructor.
     */
    explicit IIRFilter(uint32_t crossfade_samples = 0);

    /**
     * @brief   Calculate biquad cascades coefficients.
     *
     * @tparam PT   Pass type (low-pass or high-pass).
     * @param normalized_cutoff_frequency  Normalized cutoff frequency.
     */
    template<FilterPassType PT = PASS_TYPE,
            typename = std::enable_if_t<(PT == FilterPassType::LOW_PASS
                                         || PT == FilterPassType::HIGH_PASS)>>
    void calculate_cascades(double normalized_cutoff_frequency);

    /**
     * @brief   Calculate biquad cascades coefficients.
     *
     * @tparam PT   Pass type (band-pass or band-stop).
     * @param cutlow_freq  Low cutoff frequency.
     * @param cuthigh_freq  High cutoff frequency.
     */
    template<FilterPassType PT = PASS_TYPE,
            typename = std::enable_if_t<(PT == FilterPassType::BAND_PASS
                                         || PT == FilterPassType::BAND_STOP)>>
    void calculate_cascades(double cutlow_freq, double cuthigh_freq);

    CascadeFilter<PassTypeData<PASS_TYPE, N>::CASCADE_ORDER, T> _cascade_filter;

    // Pole/zero pairs without the conjugates (last one is the pair with pole on the real axis)
    PoleZeroPair _analog_pole_zero_pairs[(N + 1) / 2];

private:
    struct FrequencyConfig {
        double w1 = 0.0;
        double w2 = 0.0;
    };

    /**
     * @brief   Add a pair of S-plane pole-zero conjugates.
     *
     * @param pole_zero_pair  The pole-zero conjugates pair.
     */
    void add_pole_zero_conjugates_pair(const PoleZeroPair &pole_zero_pair);

    /**
     * @brief   Add a single S-plane pole-zero pair (pole and zero lie on the real axis).
     *
     * @param pole_zero_pair  The pole-zero pair.
     */
    void add_pole_zero_single_pair(const PoleZeroPair &pole_zero_pair);

    /**
     * @brief   Push a new block of biquad coefficients.
     *
     * @param biquad_coefficients  The biquad coefficients.
     * @return  True if pushing the coefficients was successful, false otherwise.
     */
    void push_biquad_coefficients(const BiquadCoefficients &biquad_coefficients);

    /**
     * @brief   Add a pair of S-plane pole-zero pairs.
     *
     * @note    This one is needed for band-pass/band-stop filters with an odd design order;
     *          Each real-axis pole (zero) transforms to two poles (zeros) in the complex plane after the transform,
     *          after which the poles and zeros are arranged in pairs.
     * @param pole_zero_pairs  The pole-zero pairs.
     */
    template<FilterPassType PT = PASS_TYPE,
            typename = std::enable_if_t<((PT == FilterPassType::BAND_PASS
                                          || PT == FilterPassType::BAND_STOP)
                                          && N & 1)>>
    void add_pole_zero_pairs(const std::pair<PoleZeroPair, PoleZeroPair> &pole_zero_pairs);

    // Not just CascadeFilter<N, T>::NUMBER_OF_BIQUAD_BLOCKS, because the order is doubled for band-pass/band-stop
    static constexpr uint32_t NUMBER_OF_BIQUAD_BLOCKS
            = CascadeFilter<PassTypeData<PASS_TYPE, N>::CASCADE_ORDER, T>::NUMBER_OF_BIQUAD_BLOCKS;

    PassTypeData<PASS_TYPE, N> _pass_type_data;
    FrequencyConfig _frequency_config;
};


template<uint32_t N, typename T, FilterPassType PASS_TYPE>
T IIRFilter<N, T, PASS_TYPE>::get_gain() const {
    return _cascade_filter.get_gain();
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
const T *IIRFilter<N, T, PASS_TYPE>::get_coefficients() const {
    return _cascade_filter.get_coefficients();
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
void IIRFilter<N, T, PASS_TYPE>::print_coefficients() const {
    _cascade_filter.print_coefficients();
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
T IIRFilter<N, T, PASS_TYPE>::process(const T *x, uint32_t num_samples) {
    return _cascade_filter.process(x, num_samples);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
template<typename U, typename>
void IIRFilter<N, T, PASS_TYPE>::process(const U *x, T *out, uint32_t num_samples) {
    _cascade_filter.process(x, out, num_samples);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
template<typename U, typename>
T IIRFilter<N, T, PASS_TYPE>::process(const U *x, uint32_t num_samples) {
    return _cascade_filter.process(x, num_samples);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
void IIRFilter<N, T, PASS_TYPE>::set_crossfade_samples(uint32_t crossfade_samples) {
    _cascade_filter.set_crossfade_samples(crossfade_samples);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
void IIRFilter<N, T, PASS_TYPE>::reset_state() {
    _cascade_filter.reset_state();
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
IIRFilter<N, T, PASS_TYPE>::IIRFilter(uint32_t crossfade_samples) {
    set_crossfade_samples(crossfade_samples);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType PT, typename>
void IIRFilter<N, T, PASS_TYPE>::calculate_cascades(double normalized_cutoff_frequency) {
    if (normalized_cutoff_frequency == _frequency_config.w1) {
        return;
    }

    _frequency_config.w1 = normalized_cutoff_frequency;

    _cascade_filter.set_gain(get_analog_gain());
    _cascade_filter.reset_blocks();
    _pass_type_data.init(normalized_cutoff_frequency);

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
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType PT, typename>
void IIRFilter<N, T, PASS_TYPE>::calculate_cascades(double cutlow_freq, double cuthigh_freq) {
    if (cuthigh_freq == _frequency_config.w2 && cutlow_freq == _frequency_config.w1) {
        return;
    }

    _frequency_config.w1 = cutlow_freq;
    _frequency_config.w2 = cuthigh_freq;

    _cascade_filter.set_gain(get_analog_gain());
    _cascade_filter.reset_blocks();

    constexpr double MIN_FREQ = 1e-8;
    constexpr double MAX_FREQ = 2 * M_PI - MIN_FREQ;

    if (cutlow_freq < MIN_FREQ) {
        cutlow_freq = MIN_FREQ;
    }

    if (cuthigh_freq > MAX_FREQ) {
        cuthigh_freq = MAX_FREQ;
    }

    _pass_type_data.init(cutlow_freq, cuthigh_freq);

    auto add_two_pole_zero_pairs_from_one = [&](const PoleZeroPair &pole_zero_pair) {
        auto [pole1, pole2] = _pass_type_data.transform(pole_zero_pair.pole);
        auto [zero1, zero2] = _pass_type_data.transform(pole_zero_pair.zero);

        PoleZeroPair pz1(pole1, zero1);
        add_pole_zero_conjugates_pair(pz1);

        PoleZeroPair pz2(pole2, zero2);
        add_pole_zero_conjugates_pair(pz2);
    };

    if constexpr (N & 1) {
        const PoleZeroPair &pole_zero_pair = _analog_pole_zero_pairs[(N + 1) / 2 - 1];
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
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
void IIRFilter<N, T, PASS_TYPE>::add_pole_zero_conjugates_pair(const PoleZeroPair &pole_zero_pair) {
    const Complex &pole_z = pole_zero_pair.pole;
    const Complex &zero_z = pole_zero_pair.zero;

    // (z - z0)(z - z0*)=z^2 - z(z0+z0*) + |z0|^2 = z^2 * (1 - 2*Re(z0) * z^-1 + |z0|^2 * z^-2)
    const BiquadCoefficients biquad_coefficients = {
            .b0 = 1.0,
            .b1 = -2.0 * zero_z.real(),
            .b2 = zero_z.real() * zero_z.real() + zero_z.imag() * zero_z.imag(),
            .a1 = -2.0 * pole_z.real(),
            .a2 = pole_z.real() * pole_z.real() + pole_z.imag() * pole_z.imag()
    };

    push_biquad_coefficients(biquad_coefficients);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
void IIRFilter<N, T, PASS_TYPE>::add_pole_zero_single_pair(const PoleZeroPair &pole_zero_pair) {
    const Complex &pole_z = pole_zero_pair.pole;
    const Complex &zero_z = pole_zero_pair.zero;

    // (z - z0) = z * (1 - z^-1 * z0 + 0)
    const BiquadCoefficients biquad_coefficients = {
            .b0 = 1.0,
            .b1 = -zero_z.real(),
            .b2 = 0.0,
            .a1 = -pole_z.real(),
            .a2 = 0.0
    };

    push_biquad_coefficients(biquad_coefficients);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
void IIRFilter<N, T, PASS_TYPE>::push_biquad_coefficients(const BiquadCoefficients &biquad_coefficients) {
    const double biquad_gain = _pass_type_data.calculate_gain(biquad_coefficients);
    _cascade_filter.update_biquad_gain(1.0 / biquad_gain);
    _cascade_filter.push_biquad_coefficients(biquad_coefficients);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType PT, typename>
void IIRFilter<N, T, PASS_TYPE>::add_pole_zero_pairs(const std::pair<PoleZeroPair, PoleZeroPair> &pole_zero_pairs) {
    // (z - z1)(z - z2)=z^2 - z(z1+z2) + z1z2 = z^2 * (1 - (z1+z2)z^-1 + (z1*z2)*z^-2)
    const Complex &p1 = pole_zero_pairs.first.pole;
    const Complex &p2 = pole_zero_pairs.second.pole;
    const Complex &z1 = pole_zero_pairs.first.zero;
    const Complex &z2 = pole_zero_pairs.second.zero;

    const BiquadCoefficients biquad_coefficients = {
            .b0 = 1.0,
            .b1 = -(pole_zero_pairs.first.zero.real() + pole_zero_pairs.second.zero.real()),
            .b2 = z1.imag() == 0.0
                  ? z1.real() * z2.real()
                  : std::abs(z1),
            .a1 = -(pole_zero_pairs.first.pole.real() + pole_zero_pairs.second.pole.real()),
            .a2 = p1.imag() == 0.0
                  ? p1.real() * p2.real()
                  : std::norm(p1)
    };

    push_biquad_coefficients(biquad_coefficients);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
template<typename U>
T IIRFilter<N, T, PASS_TYPE>::process(U x) {
    _cascade_filter.update_coefficients_crossfade();
    return _cascade_filter.process(x);
}

template<uint32_t N, typename T, FilterPassType PASS_TYPE>
void IIRFilter<N, T, PASS_TYPE>::process(const T *x, T *out, uint32_t num_samples) {
    _cascade_filter.update_coefficients_crossfade();
    _cascade_filter.process(x, out, num_samples);
}

}