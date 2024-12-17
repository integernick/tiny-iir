#pragma once

#include <dsp/filtering_functions.h>
#include <type_traits>

#define CASCADE_FILTER_DEBUG    1
#if CASCADE_FILTER_DEBUG > 0

#include <iostream>

#endif

namespace fast_iir {

template<unsigned int N, typename T = double>
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

    CascadeFilter() {
        update_coefficients();
    }

    explicit CascadeFilter(T *coefficients, double gain) : CascadeFilter() {
        load_biquad_coefficients(coefficients);
        set_gain(gain);
    }

    [[nodiscard]] T *get_coefficients() const {
        return _coefficients;
    }

#if CASCADE_FILTER_DEBUG > 0

    void print_coefficients() const {
        std::cout << "gain: " << _gain << std::endl;
        for (int i = 0; i < NUMBER_OF_BIQUAD_BLOCKS; ++i) {
            std::cout << "b" << i << "0: " << _coefficients[i * COEFFICIENTS_PER_BIQUAD_BLOCK]
                      << " b" << i << "1: " << _coefficients[i * COEFFICIENTS_PER_BIQUAD_BLOCK + 1]
                      << " b" << i << "2: " << _coefficients[i * COEFFICIENTS_PER_BIQUAD_BLOCK + 2]
                      << " a" << i << "0: " << _coefficients[i * COEFFICIENTS_PER_BIQUAD_BLOCK + 3]
                      << " a" << i << "1: " << _coefficients[i * COEFFICIENTS_PER_BIQUAD_BLOCK + 4]
                      << std::endl;
        }
    }

#endif

    template<typename U>
    bool push_biquad_coefficients(U b0, U b1, U b2, U a1, U a2) {
        static_assert(std::is_same<U, T>::value || std::is_floating_point<U>::value,
                      "Type mismatch: All inputs must be of the same expected type (either T or floating-point).");

        if (_num_biqaud_blocks_set >= NUMBER_OF_BIQUAD_BLOCKS) {
            return false;
        }

        U coefficients[] = {b0, b1, b2, a1, a2};

        if constexpr (std::is_same<T, U>::value) {
            memcpy(&_coefficients[_num_biqaud_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK],
                   coefficients,
                   sizeof(T) * COEFFICIENTS_PER_BIQUAD_BLOCK);
        } else if constexpr (std::is_same<T, double>::value && std::is_same<U, float>::value) {
            _coefficients[_num_biqaud_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK] = static_cast<T>(b0);
            _coefficients[_num_biqaud_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK + 1] = static_cast<T>(b1);
            _coefficients[_num_biqaud_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK + 2] = static_cast<T>(b2);
            _coefficients[_num_biqaud_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK + 3] = static_cast<T>(a1);
            _coefficients[_num_biqaud_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK + 4] = static_cast<T>(a2);
        } else if constexpr (std::is_same<T, q31_t>::value && std::is_same<U, double>::value) {
            arm_f64_to_q31(coefficients,
                           &_coefficients[_num_biqaud_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK],
                           COEFFICIENTS_PER_BIQUAD_BLOCK);
        } else if constexpr (std::is_same<T, q31_t>::value && std::is_same<U, float>::value) {
            arm_f32_to_q31(coefficients,
                           &_coefficients[_num_biqaud_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK],
                           COEFFICIENTS_PER_BIQUAD_BLOCK);
        }

        ++_num_biqaud_blocks_set;
        return true;
    }


    void update_coefficients() {
        if constexpr (std::is_same<T, float>::value) {
            arm_biquad_cascade_df2T_init_f32(&_biquad_cascade,
                                             (N + 1) / 2,
                                             _coefficients,
                                             _delay_pipeline);
        } else if constexpr (std::is_same<T, double>::value) {
            arm_biquad_cascade_df2T_init_f64(&_biquad_cascade,
                                             (N + 1) / 2,
                                             _coefficients,
                                             _delay_pipeline);
        } else if constexpr (std::is_same<T, q31_t>::value) {
            arm_biquad_cascade_df1_init_q31(&_biquad_cascade,
                                            (N + 1) / 2,
                                            _coefficients,
                                            _delay_pipeline,
                                            1);
        }
    }

    void set_gain(T gain) {
        _gain = gain;
#if CASCADE_FILTER_DEBUG > 0
        std::cout << "set_gain: " << gain << std::endl;
#endif
    }

    T process(T x) {
        x *= _gain;
        T output;

        if constexpr (std::is_same<T, double>::value) {
            arm_biquad_cascade_df2T_f64(&_biquad_cascade, &x, &output, 1);
        } else if constexpr (std::is_same<T, q31_t>::value) {
            arm_biquad_cascade_df1_q31(&_biquad_cascade, &x, &output, 1);
        }

        return output;
    }

    T process(T *x, size_t num_samples) {
        *x *= _gain;
        T output;

        if constexpr (std::is_same<T, float>::value) {
            arm_biquad_cascade_df2T_f32(&_biquad_cascade, x, _delay_pipeline, num_samples);
        } else if constexpr (std::is_same<T, double>::value) {
            arm_biquad_cascade_df2T_f64(&_biquad_cascade, x, _delay_pipeline, num_samples);
        } else if constexpr (std::is_same<T, q31_t>::value) {
            arm_biquad_cascade_df1_q31(&_biquad_cascade, x, _delay_pipeline, num_samples);
        }

        return _delay_pipeline[0];
    }

private:
    typedef typename BiquadCascade<T>::type biquad_cascade_t;

    static constexpr unsigned int NUMBER_OF_BIQUAD_BLOCKS = (N + 1) / 2;
    static constexpr unsigned int COEFFICIENTS_PER_BIQUAD_BLOCK = 5;
    static constexpr unsigned int NUMBER_OF_COEFFICIENTS
            = NUMBER_OF_BIQUAD_BLOCKS * COEFFICIENTS_PER_BIQUAD_BLOCK;

    biquad_cascade_t _biquad_cascade;

    T _gain = 0;
    T _coefficients[NUMBER_OF_COEFFICIENTS]{};
    size_t _num_biqaud_blocks_set = 0;
    T _delay_pipeline[((N + 1) / 2) * 4]{};   //TODO: 4 for direct form I
};

}
