#pragma once

#include <BiquadCascade.h>

#include <type_utils.h>

#include <iostream>

#ifndef CASCADE_FILTER_DEBUG
#define CASCADE_FILTER_DEBUG    0
#endif

namespace tiny_iir {

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
     * @brief   Initialize the filter.
     */
    void init();

    /**
     * @brief   Push a new block of biquad coefficients.
     *
     * @param biquad_coefficients  The biquad coefficients.
     * @return  True if pushing the coefficients was successful, false otherwise.
     */
    bool push_biquad_coefficients(BiquadCoefficients &biquad_coefficients);

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
    [[nodiscard]] double get_gain() const;

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
     * @param[in] x             The input buffer.
     * @param[out] out          The output buffer.
     * @param[in] num_samples   The number of samples to process.
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
     * @brief   Process a single sample (with type other than native T).
     *
     * @tparam U The type of the input sample.
     * @param x  The input sample.
     * @return   The output sample.
     */
    template<typename U, typename = std::enable_if_t<!std::is_same_v<T, U>>>
    [[nodiscard]] T process(U x);

    /**
     * @brief   Process a batch of samples (with type other than native T).
     *
     * @tparam U The type of the input samples.
     * @param[in] x             The input buffer.
     * @param[out] out          The output buffer.
     * @param[in] num_samples   The number of samples to process.
     */
    template<typename U, typename = std::enable_if_t<!std::is_same_v<T, U>>>
    void process(const U *x, T *out, uint32_t num_samples);

    /**
     * @brief   Process a batch of samples (with type other than native T).
     *
     * @tparam U The type of the input samples.
     * @param[in] x             The input buffer.
     * @param[in] num_samples   The number of samples to process.
     * @return  The final output sample.
     */
    template<typename U, typename = std::enable_if_t<!std::is_same_v<T, U>>>
    [[nodiscard]] T process(const U *x, uint32_t num_samples);

    /**
     * @brief   Update the coefficients crossfade transition.
     */
    void update_coefficients_crossfade();

#if CASCADE_FILTER_DEBUG > 0

    /**
     * @brief   Print the filter coefficients.
     */
    void print_coefficients() const;

#endif

    static constexpr uint32_t NUMBER_OF_BIQUAD_BLOCKS = (N + 1) / 2;
    static constexpr uint32_t NUMBER_OF_COEFFICIENTS
            = NUMBER_OF_BIQUAD_BLOCKS * COEFFICIENTS_PER_BIQUAD_BLOCK;
    static constexpr uint32_t DELAY_LINE_SIZE = NUMBER_OF_BIQUAD_BLOCKS * BiquadCascade<T>::BLOCK_DELAY_LINE_SIZE;

private:
    BiquadCascade<T> _biquad_cascade;

    double _gain{};
    double _gain_crossfade{};
    double _gain_delta{};

    T _coefficients[NUMBER_OF_COEFFICIENTS]{};
    uint32_t _num_biquad_blocks_set = 0;

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
        _coefficients[i * COEFFICIENTS_PER_BIQUAD_BLOCK] = 1.0;
    }

    if (_crossfade_samples > 0) {
        const double crossfade_samples_inverse = 1.0 / _crossfade_samples;
        to_native(&crossfade_samples_inverse, &_crossfade_samples_inverse, 1);
    }
}

template<uint32_t N, typename T>
const T *CascadeFilter<N, T>::get_coefficients() const {
    return _coefficients;
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::init() {
    _biquad_cascade.init(NUMBER_OF_BIQUAD_BLOCKS, _coefficients, _state);
}

template<uint32_t N, typename T>
bool CascadeFilter<N, T>::push_biquad_coefficients(BiquadCoefficients &biquad_coefficients) {
    if (_num_biquad_blocks_set >= NUMBER_OF_BIQUAD_BLOCKS) {
        return false;
    }

    T new_block[COEFFICIENTS_PER_BIQUAD_BLOCK];
    _biquad_cascade.push_biquad_coefficients(new_block, biquad_coefficients);

    const uint32_t current_block_offset = _num_biquad_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK;

    if (_crossfade_samples > 0) {
        T old_block[COEFFICIENTS_PER_BIQUAD_BLOCK];

        for (uint32_t i = 0; i < COEFFICIENTS_PER_BIQUAD_BLOCK; i++) {
            old_block[i] = _coefficients[current_block_offset + i];
        }

        for (uint32_t i = 0; i < COEFFICIENTS_PER_BIQUAD_BLOCK; i++) {
            T diff = new_block[i] - old_block[i];
            _coefficients_delta[current_block_offset + i] =
                    BiquadCascade<T>::multiply(diff, _crossfade_samples_inverse);
        }

        for (uint32_t i = 0; i < COEFFICIENTS_PER_BIQUAD_BLOCK; i++) {
            _coefficients[current_block_offset + i] = old_block[i];
        }

        _crossfade_counter = 0;
    } else {
        for (uint32_t i = 0; i < COEFFICIENTS_PER_BIQUAD_BLOCK; i++) {
            _coefficients[current_block_offset + i] = new_block[i];
        }
    }

    _num_biquad_blocks_set++;

    // TODO: fix this and remove
    // Check for NaN or infinite values in the filter state
    if (_num_biquad_blocks_set == NUMBER_OF_BIQUAD_BLOCKS) {
        for (auto &state: _state) {
            if (!std::isfinite(state)) {
                reset_state();
                break;
            }
        }
    }

    return true;
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::update_biquad_gain(double biquad_gain_inv) {
    std::cout << "biquad_gain_inv: " << biquad_gain_inv << std::endl;

    _gain *= biquad_gain_inv;

    std::cout << "_gain: " << _gain << std::endl;

    if (_crossfade_samples > 0) {
        _gain_delta = _gain - _gain_crossfade;
        _gain_delta = BiquadCascade<T>::multiply(_gain_delta, _crossfade_samples_inverse);
    }
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::reset_state() {
    memset(_state, 0, sizeof(_state));
    _is_set_once = false;
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::reset_blocks() {
    _num_biquad_blocks_set = 0;
}

template<uint32_t N, typename T>
double CascadeFilter<N, T>::get_gain() const {
    return _gain;
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::set_gain(double gain) {
    _gain = gain;
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::set_crossfade_samples(uint32_t crossfade_samples) {
    _crossfade_samples = crossfade_samples;
}

template<uint32_t N, typename T>
T CascadeFilter<N, T>::process(T x) {
    double x_double;
    to_double(&x, &x_double, 1);
    x_double *= _gain_crossfade;
    to_native(&x_double, &x, 1);

    T output = 0;

    _biquad_cascade.process_cascade(&x, &output, 1);

    return output;
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::process(const T *x, T *out, uint32_t num_samples) {
    if (num_samples == 0) {
        return;
    }

    T input_buffer[num_samples];

    if constexpr (std::is_same_v<T, double> or std::is_same_v<T, float>) {
        for (uint32_t i = 0; i < num_samples; ++i) {
            input_buffer[i] = x[i] * _gain_crossfade;
        }
    } else {
        for (uint32_t i = 0; i < num_samples; ++i) {
            double x_double;
            to_double(&x[i], &x_double, 1);
            x_double *= _gain_crossfade;
            to_native(&x_double, &input_buffer[i], 1);
        }
    }

    for (uint32_t i = 0; i < num_samples; ++i) {
        input_buffer[i] = BiquadCascade<T>::multiply(x[i], _gain_crossfade);
    }

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
T CascadeFilter<N, T>::process(U x) {
    T x_native;

    if constexpr (std::is_same_v<U, double> or std::is_same_v<U, float>) {
        x *= _gain_crossfade;
        to_native(&x, &x_native, 1);
    } else {
        double x_double;
        to_double(&x, &x_double, 1);
        x_double *= _gain_crossfade;
        to_native(&x_double, &x_native, 1);
    }

    T output = 0;

    _biquad_cascade.process_cascade(&x_native, &output, 1);

    return output;
}

template<uint32_t N, typename T>
template<typename U, typename>
void CascadeFilter<N, T>::process(const U *x, T *out, uint32_t num_samples) {
    if (num_samples == 0) {
        return;
    }

    T x_native[num_samples];

    if constexpr (std::is_same_v<U, double> or std::is_same_v<U, float>) {
        for (uint32_t i = 0; i < num_samples; ++i) {
            const double sample = x[i] * _gain_crossfade;
            to_native(&sample, &x_native[i], 1);
        }
    } else {
        for (uint32_t i = 0; i < num_samples; ++i) {
            double x_double;
            to_double(&x[i], &x_double, 1);
            x_double *= _gain_crossfade;
            to_native(&x_double, &x_native[i], 1);
        }
    }

    _biquad_cascade.process_cascade(x_native, out, num_samples);
}

template<uint32_t N, typename T>
template<typename U, typename>
T CascadeFilter<N, T>::process(const U *x, uint32_t num_samples) {
    if (num_samples == 0) {
        return 0;
    }

    T output_buffer[num_samples];
    process(x, output_buffer, num_samples);

    return output_buffer[num_samples - 1];
}

template<uint32_t N, typename T>
void CascadeFilter<N, T>::update_coefficients_crossfade() {
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

        _crossfade_counter += 1;
    }
}

#if CASCADE_FILTER_DEBUG > 0
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
        double biquad_coefficients_double[COEFFICIENTS_PER_BIQUAD_BLOCK];
        to_double(&_coefficients[COEFFICIENTS_PER_BIQUAD_BLOCK * i],
                  biquad_coefficients_double, COEFFICIENTS_PER_BIQUAD_BLOCK);
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
