#pragma once

#include "utils.h"
#include "BiquadCascade.h"

#define CASCADE_FILTER_DEBUG    0

namespace tiny_iir {

template<uint32_t N, typename T, uint32_t CROSSFADE_SAMPLES>
class CascadeFilter {
public:
    CascadeFilter();

    [[nodiscard]] const T *get_coefficients() const;

    bool push_biquad_coefficients(double b0, double b1, double b2, double a1, double a2);

    void reset_state();

    void reset();

    [[nodiscard]] T get_gain() const;

    void set_gain(double gain);

    [[nodiscard]] T process(T x);

    void process(const T *x, T *out, uint32_t num_samples);

    [[nodiscard]] T process(const T *x, uint32_t num_samples);

    template<typename U, typename std::enable_if_t<!std::is_same_v<T, U>, bool> = 0>
    [[nodiscard]] T process(U x);

    template<typename U, typename std::enable_if_t<!std::is_same_v<T, U>, bool> = 0>
    void process(const U *x, T *out, uint32_t num_samples);

    template<typename U, typename std::enable_if_t<!std::is_same_v<T, U>, bool> = 0>
    [[nodiscard]] T process(const U *x, uint32_t num_samples);

#if CASCADE_FILTER_DEBUG > 0

    void print_coefficients() const;

#endif

    static constexpr uint32_t NUMBER_OF_BIQUAD_BLOCKS = (N + 1) / 2;
    static constexpr uint32_t NUMBER_OF_COEFFICIENTS
            = NUMBER_OF_BIQUAD_BLOCKS * COEFFICIENTS_PER_BIQUAD_BLOCK;
    static constexpr uint32_t DELAY_LINE_SIZE = std::is_same<T, q31_t>::value
                                                    ? NUMBER_OF_BIQUAD_BLOCKS * 4
                                                    : NUMBER_OF_BIQUAD_BLOCKS * 2;

private:
    typedef typename BiquadCascade<T>::type biquad_cascade_t;

    template<typename = std::enable_if<(false)>>
    void update_coefficients_crossfade();

    biquad_cascade_t _biquad_cascade;

    T _gain{};
    T _coefficients[NUMBER_OF_COEFFICIENTS]{};
    uint32_t _num_biquad_blocks_set = 0;

    T _state[DELAY_LINE_SIZE]{};

    // Smooth transition (crossfade) between old and new coefficients
    T _coefficients_delta[NUMBER_OF_COEFFICIENTS]{};
    uint32_t _crossfaded_samples = 0;
    T _crossfade_samples_inverse{};
};

template<uint32_t N, typename T, uint32_t CROSSFADE_SAMPLES>
CascadeFilter<N, T, CROSSFADE_SAMPLES>::CascadeFilter() {
    // Initialize as passthrough filter
    for (uint32_t i = 0; i < NUMBER_OF_BIQUAD_BLOCKS; ++i) {
        _coefficients[i * COEFFICIENTS_PER_BIQUAD_BLOCK] = 1.0;
    }

    BiquadCascade<T>::init(&_biquad_cascade, NUMBER_OF_BIQUAD_BLOCKS, _coefficients, _state);

    if constexpr (CROSSFADE_SAMPLES > 0) {
        const double crossfade_samples_inverse = 1.0 / CROSSFADE_SAMPLES;
        to_native(&crossfade_samples_inverse, &_crossfade_samples_inverse, 1);
    }
}

template<uint32_t N, typename T, uint32_t CROSSFADE_SAMPLES>
const T *CascadeFilter<N, T, CROSSFADE_SAMPLES>::get_coefficients() const {
    return _coefficients;
}

template<uint32_t N, typename T, uint32_t CROSSFADE_SAMPLES>
bool CascadeFilter<N, T, CROSSFADE_SAMPLES>::push_biquad_coefficients(double b0, double b1, double b2,
                                                                      double a1, double a2) {
    if (_num_biquad_blocks_set >= NUMBER_OF_BIQUAD_BLOCKS) {
        return false;
    }

    T new_block[COEFFICIENTS_PER_BIQUAD_BLOCK];
    BiquadCascade<T>::push_biquad_coefficients(new_block, b0, b1, b2, a1, a2);

    const uint32_t current_block_offset = _num_biquad_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK;

    if constexpr (CROSSFADE_SAMPLES > 0) {
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

        _crossfaded_samples = 0;
    } else {
        for (uint32_t i = 0; i < COEFFICIENTS_PER_BIQUAD_BLOCK; i++) {
            _coefficients[current_block_offset + i] = new_block[i];
        }
    }

    _num_biquad_blocks_set++;

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

template<uint32_t N, typename T, uint32_t CROSSFADE_SAMPLES>
void CascadeFilter<N, T, CROSSFADE_SAMPLES>::reset_state() {
    memset(_state, 0, sizeof(_state));
}

template<uint32_t N, typename T, uint32_t CROSSFADE_SAMPLES>
void CascadeFilter<N, T, CROSSFADE_SAMPLES>::reset() {
    _num_biquad_blocks_set = 0;
}

template<uint32_t N, typename T, uint32_t CROSSFADE_SAMPLES>
T CascadeFilter<N, T, CROSSFADE_SAMPLES>::get_gain() const {
    return _gain;
}

template<uint32_t N, typename T, uint32_t CROSSFADE_SAMPLES>
void CascadeFilter<N, T, CROSSFADE_SAMPLES>::set_gain(double gain) {
    if constexpr (std::is_same<T, double>::value) {
        _gain = gain;
    } else {
        to_native(&gain, &_gain, 1);
    }
}

template<uint32_t N, typename T, uint32_t CROSSFADE_SAMPLES>
T CascadeFilter<N, T, CROSSFADE_SAMPLES>::process(T x) {
    x = BiquadCascade<T>::multiply(x, _gain);
    T output = 0;

    if constexpr (CROSSFADE_SAMPLES > 0) {
        update_coefficients_crossfade();
    }

    BiquadCascade<T>::process_cascade(&_biquad_cascade, &x, &output, 1);

    return output;
}

template<uint32_t N, typename T, uint32_t CROSSFADE_SAMPLES>
void CascadeFilter<N, T, CROSSFADE_SAMPLES>::process(const T *x, T *out, uint32_t num_samples) {
    if (num_samples == 0) {
        return;
    }
    T input_buffer[num_samples];
    for (uint32_t i = 0; i < num_samples; ++i) {
        input_buffer[i] = BiquadCascade<T>::multiply(x[i], _gain);
    }

    if constexpr (CROSSFADE_SAMPLES > 0) {
        update_coefficients_crossfade();
    }

    BiquadCascade<T>::process_cascade(&_biquad_cascade, input_buffer, out, num_samples);
}

template<uint32_t N, typename T, uint32_t CROSSFADE_SAMPLES>
T CascadeFilter<N, T, CROSSFADE_SAMPLES>::process(const T *x, uint32_t num_samples) {
    if (num_samples == 0) {
        return 0;
    }

    T output_buffer[num_samples];
    process(x, output_buffer, num_samples);

    return output_buffer[num_samples - 1];
}

template<uint32_t N, typename T, uint32_t CROSSFADE_SAMPLES>
template<typename U, std::enable_if_t<!std::is_same_v<T, U>, bool>>
T CascadeFilter<N, T, CROSSFADE_SAMPLES>::process(U x) {
    T x_native;
    to_native(&x, &x_native, 1);
    x_native = BiquadCascade<T>::multiply(x_native, _gain);
    T output = 0;

    if constexpr (CROSSFADE_SAMPLES > 0) {
        update_coefficients_crossfade();
    }

    BiquadCascade<T>::process_cascade(&_biquad_cascade, &x_native, &output, 1);
    return output;
}

template<uint32_t N, typename T, uint32_t CROSSFADE_SAMPLES>
template<typename U, std::enable_if_t<!std::is_same_v<T, U>, bool>>
void CascadeFilter<N, T, CROSSFADE_SAMPLES>::process(const U *x, T *out, uint32_t num_samples) {
    if (num_samples == 0) {
        return;
    }

    T x_native[num_samples];
    to_native(x, x_native, num_samples);

    for (uint32_t i = 0; i < num_samples; ++i) {
        x_native[i] = BiquadCascade<T>::multiply(x_native[i], _gain);
    }

    if constexpr (CROSSFADE_SAMPLES > 0) {
        update_coefficients_crossfade();
    }

    BiquadCascade<T>::process_cascade(&_biquad_cascade, x_native, out, num_samples);
}

template<uint32_t N, typename T, uint32_t CROSSFADE_SAMPLES>
template<typename U, std::enable_if_t<!std::is_same_v<T, U>, bool>>
T CascadeFilter<N, T, CROSSFADE_SAMPLES>::process(const U *x, uint32_t num_samples) {
    if (num_samples == 0) {
        return 0;
    }

    T output_buffer[num_samples];
    process(x, output_buffer, num_samples);

    return output_buffer[num_samples - 1];
}

#if CASCADE_FILTER_DEBUG > 0
template<uint32_t N, typename T, uint32_t CROSSFADE_SAMPLES>
void CascadeFilter<N, T, CROSSFADE_SAMPLES>::print_coefficients() const {
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
                 "1 %.15f %.15f 1 %.15f %.15f\n",
                 biquad_coefficients_double[1],
                 biquad_coefficients_double[2],
                 -biquad_coefficients_double[3],
                 -biquad_coefficients_double[4]);
        printf("%s", debug_buf);
    }
}
#endif

template<uint32_t N, typename T, uint32_t CROSSFADE_SAMPLES>
template<typename>
void CascadeFilter<N, T, CROSSFADE_SAMPLES>::update_coefficients_crossfade() {
    if (_crossfaded_samples < CROSSFADE_SAMPLES) {
        for (uint32_t i = 0; i < NUMBER_OF_COEFFICIENTS; ++i) {
            _coefficients[i] += _coefficients_delta[i];
        }

        _crossfaded_samples += 1;
    }
}

} // namespace tiny_iir
