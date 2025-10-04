#pragma once

#include <BiquadCascade.h>

#include <type_utils.h>

#include <cstring>

#ifndef TINY_IIR_CASCADE_FILTER_DEBUG
#define TINY_IIR_CASCADE_FILTER_DEBUG    0
#endif

static_assert(TINY_IIR_CHUNK_SIZE > 0, "TINY_IIR_CHUNK_SIZE must be > 0");

namespace tiny_iir {

/**
 * @brief   Cascade filter with crossfade functionality.
 *
 * @details Implements the transfer function of series biquad blocks
 *          and optional crossfade between old and new coefficients.
 * @tparam N   Analog prototype filter order.
 * @tparam T   Data type.
 */
template<uint32_t N, typename T>
class CascadeFilter {
public:
    /**
     * @brief   Constructor.
     */
    CascadeFilter();

    /**
     * @brief   Get the filter coefficients.
     *
     * @return  Pointer to the coefficients array.
     */
    [[nodiscard]] const T *get_coefficients() const;

    /**
     * @brief   Initialize the biquad cascade coefficients.
     */
    void init_coefficients();

    /**
     * @brief   Initialize the filter.
     *
     * @details Clears the state buffer.
     */
    void init();

    /**
     * @brief   Reset the cascade state.
     */
    void reset();

    /**
     * @brief   Push a new block of biquad coefficients.
     *
     * @param biquad_coefficients  The biquad coefficients.
     * @return  True if pushing the coefficients was successful, false otherwise.
     */
    bool push_biquad_coefficients(BiquadCoefficients &biquad_coefficients);

    /**
     * @brief   Get the biquad coefficients (double representation).
     *
     * @param biquad_idx  The index of the biquad.
     * @return  The biquad coefficients.
     */
    [[nodiscard]] BiquadCoefficients get_biquad_coefficients(uint32_t biquad_idx) const;

    /**
     * @brief   Update total filter gain with a new biquad gain.
     *
     * @note    This method is used by the external IIRFilter class that calculates gain depending on the pass type.
     * @param biquad_gain_inv  The new biquad gain inverse.
     */
    void update_biquad_gain(double biquad_gain_inv);

    /**
     * @brief   Reset the filter state.
     */
    void reset_state();

    /**
     * @brief   Reset the blocks number.
     */
    void reset_blocks();

    /**
     * @brief   Get the filter gain.
     *
     * @return  The filter gain.
     */
    [[nodiscard]] T get_gain() const;

    /**
     * @brief   Set the filter gain.
     *
     * @param gain  The filter gain.
     */
    void set_gain(double gain);

    /**
     * @brief   Set the number of samples to smooth the transition between old and new coefficients.
     *
     * @param crossfade_samples  The number of samples to smooth the transition between old and new coefficients.
     */
    void set_crossfade_samples(uint32_t crossfade_samples);

    /**
     * @brief   Process a single sample.
     *
     * @param x  The input sample.
     * @return   The output sample.
     */
    [[nodiscard]] T process(T x);

    /**
     * @brief   Process a batch of samples.
     *
     * @param[in] in            The input buffer.
     * @param[out] out          The output buffer.
     * @param[in] num_samples   The number of samples to process.
     */
    void process(const T *in, T *out, uint32_t num_samples);

    /**
     * @brief   Process a batch of samples.
     *
     * @param[in] x         The input buffer.
     * @param num_samples   The number of samples to process.
     * @return  The final output sample.
     */
    [[nodiscard]] T process(const T *x, uint32_t num_samples);

    /**
     * @brief   Process a single sample (with type other than native T).
     *
     * @tparam U The type of the input sample.
     * @param x  The input sample.
     * @return   The output sample.
     */
    template<typename U, typename = std::enable_if_t<!std::is_same_v<T, U>>>
    [[nodiscard]] U process(U x);

    /**
     * @brief   Process a batch of samples (with type other than native T).
     *
     * @tparam U The type of the input samples.
     * @param[in] in            The input buffer.
     * @param[out] out          The output buffer.
     * @param[in] num_samples   The number of samples to process.
     */
    template<typename U, typename = std::enable_if_t<!std::is_same_v<T, U>>>
    void process(const U *in, U *out, uint32_t num_samples);

    /**
     * @brief   Process a batch of samples (with type other than native T).
     *
     * @tparam U The type of the input samples.
     * @param[in] in            The input buffer.
     * @param[in] num_samples   The number of samples to process.
     * @return  The final output sample.
     */
    template<typename U, typename = std::enable_if_t<!std::is_same_v<T, U>>>
    [[nodiscard]] U process(const U *in, uint32_t num_samples);

    /**
     * @brief   Update the coefficients crossfade transition.
     */
    void update_coefficients_crossfade();

#if TINY_IIR_CASCADE_FILTER_DEBUG > 0

    /**
     * @brief   Print the filter coefficients.
     */
    void print_coefficients() const;

#endif

    static constexpr uint32_t NUMBER_OF_BIQUAD_BLOCKS = (N + 1) / 2;
    static constexpr uint32_t NUMBER_OF_COEFFICIENTS
            = NUMBER_OF_BIQUAD_BLOCKS * coeffs_per_stage<T>::value;
    static constexpr uint32_t DELAY_LINE_SIZE = NUMBER_OF_BIQUAD_BLOCKS * BiquadCascade<T>::BLOCK_DELAY_LINE_SIZE;
    static constexpr uint32_t PROCESS_CHUNK_SIZE = TINY_IIR_CHUNK_SIZE;

private:
    BiquadCascade<T> _biquad_cascade;

    T _gain{BiquadCascade<T>::UNITY};
    T _gain_crossfade{BiquadCascade<T>::UNITY};
    T _gain_delta{};

    T _coefficients[NUMBER_OF_COEFFICIENTS]{};

    T _state[DELAY_LINE_SIZE]{};

    // Smoothen transition between old and new coefficients (crossfade)
    uint32_t _crossfade_samples = 0;
    uint32_t _crossfade_counter = 0;
    T _coefficients_delta[NUMBER_OF_COEFFICIENTS]{};
    T _crossfade_samples_inverse{};

    bool _is_set_once = false;
};

template<uint32_t N, typename T>
CascadeFilter<N, T>::CascadeFilter() {
    // Initialize as passthrough filter
    for (uint32_t i = 0; i < NUMBER_OF_BIQUAD_BLOCKS; ++i) {
        _coefficients[i * coeffs_per_stage<T>::value] = 1.0;
    }
}

template<uint32_t N, typename T>
const T *CascadeFilter<N, T>::get_coefficients() const {
    return _coefficients;
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::init_coefficients() {
    _biquad_cascade.set_coefficients(_coefficients);
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::init() {
    _biquad_cascade.set_coefficients(_coefficients);

    if constexpr (std::is_same_v<T, q31_t> || std::is_same_v<T, q15_t>) {
        _biquad_cascade.init(NUMBER_OF_BIQUAD_BLOCKS, _state);
    } else {
        if (not _is_set_once) {
            _biquad_cascade.init(NUMBER_OF_BIQUAD_BLOCKS, _state);
        }
    }
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::reset() {
    // Only reset internal counters/state of the biquad cascade; do not touch coefficients
    _biquad_cascade.reset();
}

template<uint32_t N, typename T>
bool CascadeFilter<N, T>::push_biquad_coefficients(BiquadCoefficients &biquad_coefficients) {
    if (_biquad_cascade.get_num_biquad_blocks_set() >= NUMBER_OF_BIQUAD_BLOCKS) {
        return false;
    }

    T new_block[coeffs_per_stage<T>::value];
    _biquad_cascade.push_biquad_coefficients(new_block, biquad_coefficients);

    const uint32_t current_block_offset
            = (_biquad_cascade.get_num_biquad_blocks_set() - 1) * coeffs_per_stage<T>::value;

    if (_crossfade_samples > 0) {
        for (uint32_t i = 0; i < coeffs_per_stage<T>::value; ++i) {
            const uint32_t idx = current_block_offset + i;
            const T diff = new_block[i] - _coefficients[idx];
            scale(&diff, _crossfade_samples_inverse, &_coefficients_delta[idx], 1);
        }

        _crossfade_counter = 0;
    } else {
        memcpy(_coefficients + current_block_offset, new_block, sizeof(T) * coeffs_per_stage<T>::value);
    }

    return true;
}

template<uint32_t N, typename T>
BiquadCoefficients CascadeFilter<N, T>::get_biquad_coefficients(uint32_t biquad_idx) const {
    return _biquad_cascade.get_biquad_coefficients(biquad_idx);
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::update_biquad_gain(double biquad_gain_inv) {
    _biquad_cascade.update_gain(biquad_gain_inv, _gain);

    if (_crossfade_samples > 0) {
        _gain_delta = _gain - _gain_crossfade;
        scale(&_gain_delta, _crossfade_samples_inverse, &_gain_delta, 1);
    } else {
        _gain_crossfade = _gain;
    }
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::reset_state() {
    memset(_state, 0, sizeof(_state));
    _is_set_once = false;
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::reset_blocks() {
    _biquad_cascade.reset();
}

template<uint32_t N, typename T>
T CascadeFilter<N, T>::get_gain() const {
    return _gain;
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::set_gain(double gain) {
    to_native(&gain, &_gain, 1);
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::set_crossfade_samples(uint32_t crossfade_samples) {
    _crossfade_samples = crossfade_samples;

    if (_crossfade_samples > 0) {
        double inv = 1.0 / static_cast<double>(_crossfade_samples);
        to_native(&inv, &_crossfade_samples_inverse, 1);
    } else {
        _crossfade_samples_inverse = {};
    }
}

template<uint32_t N, typename T>
T CascadeFilter<N, T>::process(T x) {
    scale(&x, _gain_crossfade, &x, 1);
    T output = 0;
    _biquad_cascade.process_cascade(&x, &output, 1);

    return output;
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::process(const T *in, T *out, uint32_t num_samples) {
    if (num_samples == 0) {
        return;
    }

    T input_buffer[num_samples];
    scale(in, _gain_crossfade, input_buffer, num_samples);

    _biquad_cascade.process_cascade(input_buffer, out, num_samples);
}

template<uint32_t N, typename T>
T CascadeFilter<N, T>::process(const T *x, uint32_t num_samples) {
    if (num_samples == 0) {
        return 0;
    }

    T output_buffer[num_samples];
    process(x, output_buffer, num_samples);

    return output_buffer[num_samples - 1];
}

template<uint32_t N, typename T>
template<typename U, typename>
U CascadeFilter<N, T>::process(U x) {
    T x_native;
    to_native(&x, &x_native, 1);

    T out_native = process(x_native);
    U out_derived;
    to_native<U, T>(&out_native, &out_derived, 1);
    return out_derived;
}

template<uint32_t N, typename T>
template<typename U, typename>
void CascadeFilter<N, T>::process(const U *in, U *out, uint32_t num_samples) {
    if (num_samples == 0) {
        return;
    }

    T in_native[PROCESS_CHUNK_SIZE];

    for (uint32_t i = 0; i < num_samples; ) {
        const uint32_t n = std::min<uint32_t>(PROCESS_CHUNK_SIZE, num_samples - i);
        to_native(in + i, in_native, n);
        process(in_native, in_native, n);
        to_native<U, T>(in_native, out + i, n);
        i += n;
    }
}

template<uint32_t N, typename T>
template<typename U, typename>
U CascadeFilter<N, T>::process(const U *in, uint32_t num_samples) {
    if (num_samples == 0) {
        return 0;
    }

    U last{};
    T in_native[PROCESS_CHUNK_SIZE];

    for (uint32_t i = 0; i < num_samples; ) {
        const uint32_t n = std::min<uint32_t>(PROCESS_CHUNK_SIZE, num_samples - i);
        to_native(in + i, in_native, n);
        process(in_native, in_native, n);
        to_native<U, T>(&in_native[n - 1], &last, 1);
        i += n;
    }

    return last;
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::update_coefficients_crossfade() {
    if (_crossfade_samples == 0) {
        _gain_crossfade = _gain;
        return;
    }

    if (not _is_set_once) {
        _gain_crossfade = _gain;

        for (uint32_t i = 0; i < NUMBER_OF_COEFFICIENTS; ++i) {
            _coefficients[i] += _coefficients_delta[i] * _crossfade_samples;
        }

        _crossfade_counter = _crossfade_samples;
        _is_set_once = true;
        return;
    }

    if (_crossfade_counter < _crossfade_samples) {
        _gain_crossfade += _gain_delta;

        for (uint32_t i = 0; i < NUMBER_OF_COEFFICIENTS; ++i) {
            _coefficients[i] += _coefficients_delta[i];
        }

        _crossfade_counter++;
    }
}

#if TINY_IIR_CASCADE_FILTER_DEBUG > 0

template<uint32_t N, typename T>
void CascadeFilter<N, T>::print_coefficients() const {
    double gain_double;
    to_double(&_gain, &gain_double, 1);
    char debug_buf[128];
    snprintf(debug_buf, sizeof(debug_buf), "gain: %.15f\n", gain_double);
    printf("%s", debug_buf);
    snprintf(debug_buf, sizeof(debug_buf), "coefficients:\n");
    printf("%s", debug_buf);
    for (int i = 0; i < NUMBER_OF_BIQUAD_BLOCKS; ++i) {
        double biquad_coefficients_double[coeffs_per_stage<T>::value];
        to_double(&_coefficients[coeffs_per_stage<T>::value * i],
                  biquad_coefficients_double, coeffs_per_stage<T>::value);
        snprintf(debug_buf, sizeof(debug_buf),
                 "%.15f %.15f %.15f 1 %.15f %.15f\n",
                 biquad_coefficients_double[0],
                 biquad_coefficients_double[1],
                 biquad_coefficients_double[2],
                 -biquad_coefficients_double[3],
                 -biquad_coefficients_double[4]);
        printf("%s", debug_buf);
    }
}

#endif

} // namespace tiny_iir
