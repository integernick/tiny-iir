#pragma once

#include <dsp/filtering_functions.h>
#include <type_traits>

#define CASCADE_FILTER_DEBUG    1
#if CASCADE_FILTER_DEBUG > 0
#include <iostream>
#include <iomanip>
#endif

namespace tiny_iir {

template<unsigned int N, typename T = double, unsigned int MAX_SAMPLES = 8>
class CascadeFilter {
public:
    template<typename>
    struct BiquadCascade;

    template<>
    struct BiquadCascade<float> {
        typedef arm_biquad_cascade_df2T_instance_f32 type;
    };

    template<>
    struct BiquadCascade<double> {
        typedef arm_biquad_cascade_df2T_instance_f64 type;
    };

    template<>
    struct BiquadCascade<q31_t> {
        typedef arm_biquad_casd_df1_inst_q31 type;
    };

    CascadeFilter() = default;

    explicit CascadeFilter(T *coefficients, double gain) : CascadeFilter() {
        load_biquad_coefficients(coefficients);
        set_gain(gain);
    }

    [[nodiscard]] const T *get_coefficients() const {
        return _coefficients;
    }

#if CASCADE_FILTER_DEBUG > 0
    void print_coefficients() const {
        std::cout << std::setprecision(15) << "gain: " << _gain << std::endl;
        for (int i = 0; i < NUMBER_OF_BIQUAD_BLOCKS; ++i) {
            std::cout << std::setprecision(15)
                      << _coefficients[i * COEFFICIENTS_PER_BIQUAD_BLOCK]
                      << ", " << _coefficients[i * COEFFICIENTS_PER_BIQUAD_BLOCK + 1]
                      << ", " << _coefficients[i * COEFFICIENTS_PER_BIQUAD_BLOCK + 2]
                      << ", 1"
                      << ", " << _coefficients[i * COEFFICIENTS_PER_BIQUAD_BLOCK + 3]
                      << ", " << _coefficients[i * COEFFICIENTS_PER_BIQUAD_BLOCK + 4] << ";"
                      << std::endl;
        }
    }
#endif

    template<typename U>
    bool push_biquad_coefficients(U b0, U b1, U b2, U a1, U a2) {
        static_assert(std::is_same<U, T>::value || std::is_floating_point<U>::value,
                      "Type mismatch: All inputs must be of the same expected type (either T or floating-point).");

        if (_num_biquad_blocks_set >= NUMBER_OF_BIQUAD_BLOCKS) {
            return false;
        }

        U coefficients[] = {b0, b1, b2, -a1, -a2};

        if constexpr (std::is_same<T, U>::value) {
            memcpy(&_coefficients[_num_biquad_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK],
                   coefficients,
                   sizeof(T) * COEFFICIENTS_PER_BIQUAD_BLOCK);
        } else if constexpr (std::is_same<T, double>::value && std::is_same<U, float>::value) {
            _coefficients[_num_biquad_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK] = static_cast<T>(b0);
            _coefficients[_num_biquad_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK + 1] = static_cast<T>(b1);
            _coefficients[_num_biquad_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK + 2] = static_cast<T>(b2);
            _coefficients[_num_biquad_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK + 3] = static_cast<T>(-a1);
            _coefficients[_num_biquad_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK + 4] = static_cast<T>(-a2);
        } else if constexpr (std::is_same<T, q31_t>::value && std::is_same<U, double>::value) {
            arm_f64_to_q31(coefficients,
                           &_coefficients[_num_biquad_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK],
                           COEFFICIENTS_PER_BIQUAD_BLOCK);
        } else if constexpr (std::is_same<T, q31_t>::value && std::is_same<U, float>::value) {
            arm_f32_to_q31(coefficients,
                           &_coefficients[_num_biquad_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK],
                           COEFFICIENTS_PER_BIQUAD_BLOCK);
        }

        ++_num_biquad_blocks_set;
        return true;
    }

    void init_biquad_cascades() {
        if constexpr (std::is_same<T, float>::value) {
            arm_biquad_cascade_df2T_init_f32(&_biquad_cascade,
                                             _num_biquad_blocks_set,
                                             _coefficients,
                                             _delay_pipeline);
        } else if constexpr (std::is_same<T, double>::value) {
            arm_biquad_cascade_df2T_init_f64(&_biquad_cascade,
                                             _num_biquad_blocks_set,
                                             _coefficients,
                                             _delay_pipeline);
        } else if constexpr (std::is_same<T, q31_t>::value) {
            arm_biquad_cascade_df1_init_q31(&_biquad_cascade,
                                            _num_biquad_blocks_set,
                                            _coefficients,
                                            _delay_pipeline,
                                            1);
        }
    }

    void reset() {
        _num_biquad_blocks_set = 0;
        memset(_coefficients, 0, sizeof(_coefficients));
        memset(_delay_pipeline, 0, sizeof(_delay_pipeline));
        memset(_output_buffer, 0, sizeof(_output_buffer));
    }

    [[nodiscard]] T get_gain() const {
        return _gain;
    }

    void set_gain(T gain) {
        _gain = gain;
    }

    T process(T x) {
        x *= _gain;
        T output = 0;
        arm_biquad_cascade(&_biquad_cascade, &x, &output, 1);
        return output;
    }

    T process(T *x, size_t num_samples) {
        static_assert(num_samples <= MAX_SAMPLES, "Output buffer size too small");
        for (size_t i = 0; i < num_samples; ++i) {
            x[i] *= _gain;
        }
        arm_biquad_cascade(&_biquad_cascade, x, _output_buffer, num_samples);
        return _output_buffer[num_samples - 1];
    }

private:
    typedef typename BiquadCascade<T>::type biquad_cascade_t;

    void arm_biquad_cascade(biquad_cascade_t *biquad_cascade, const T *src, T *dst, uint32_t block_size) {
        if constexpr (std::is_same<T, float>::value) {
            arm_biquad_cascade_df2T_f32(biquad_cascade, src, dst, block_size);
        } else if constexpr (std::is_same<T, double>::value) {
            arm_biquad_cascade_df2T_f64(biquad_cascade, src, dst, block_size);
        } else if constexpr (std::is_same<T, q31_t>::value) {
            arm_biquad_cascade_df1_q31(biquad_cascade, src, dst, block_size);
        }
    }

    static constexpr unsigned int NUMBER_OF_BIQUAD_BLOCKS = (N + 1) / 2;
    static constexpr unsigned int COEFFICIENTS_PER_BIQUAD_BLOCK = 5;
    static constexpr unsigned int NUMBER_OF_COEFFICIENTS
            = NUMBER_OF_BIQUAD_BLOCKS * COEFFICIENTS_PER_BIQUAD_BLOCK;
    static constexpr unsigned int DELAY_LINE_SIZE = std::is_same<T, q31_t>::value
            ? NUMBER_OF_BIQUAD_BLOCKS * 4
            : NUMBER_OF_BIQUAD_BLOCKS * 2;

    biquad_cascade_t _biquad_cascade;

    T _gain = 1.0;
    T _coefficients[NUMBER_OF_COEFFICIENTS]{};
    size_t _num_biquad_blocks_set = 0;

    T _delay_pipeline[DELAY_LINE_SIZE]{};

    T _output_buffer[MAX_SAMPLES]{};
};

}
