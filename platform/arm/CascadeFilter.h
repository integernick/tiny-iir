#pragma once

#include "utils.h"

#include <dsp/filtering_functions.h>
#include <type_traits>

#define CASCADE_FILTER_DEBUG    0
#if CASCADE_FILTER_DEBUG > 0

#include <cstdio>

#endif

namespace tiny_iir {

template<unsigned int N, typename T = double>
class CascadeFilter {
public:
    template<typename>
    struct BiquadCascade;

    template<>
    struct BiquadCascade<float> {
        typedef arm_biquad_cascade_df2T_instance_f32 type;

        static void init(type *cascade, uint8_t num_stages, float *coefficients, float *delay_pipeline) {
            arm_biquad_cascade_df2T_init_f32(cascade, num_stages, coefficients, delay_pipeline);
        }

        static void process_cascade(type *cascade, float *src, float *dst, uint32_t block_size) {
            arm_biquad_cascade_df2T_f32(cascade, src, dst, block_size);
        }

        [[nodiscard]] static float multiply(float x, float y) {
            return x * y;
        }

        static void
        push_biquad_coefficients(float *coefficients, double b0, double b1, double b2, double a1, double a2) {
            coefficients[0] = static_cast<float>(b0);
            coefficients[1] = static_cast<float>(b1);
            coefficients[2] = static_cast<float>(b2);
            coefficients[3] = static_cast<float>(-a1);
            coefficients[4] = static_cast<float>(-a2);
        }
    };

    template<>
    struct BiquadCascade<double> {
        typedef arm_biquad_cascade_df2T_instance_f64 type;

        static void init(type *cascade, uint8_t num_stages, double *coefficients, double *delay_pipeline) {
            arm_biquad_cascade_df2T_init_f64(cascade, num_stages, coefficients, delay_pipeline);
        }

        static void process_cascade(type *cascade, double *src, double *dst, uint32_t block_size) {
            arm_biquad_cascade_df2T_f64(cascade, src, dst, block_size);
        }

        [[nodiscard]] static double multiply(double x, double y) {
            return x * y;
        }

        static void
        push_biquad_coefficients(double *coefficients, double b0, double b1, double b2, double a1, double a2) {
            coefficients[0] = b0;
            coefficients[1] = b1;
            coefficients[2] = b2;
            coefficients[3] = -a1;
            coefficients[4] = -a2;
        }
    };

    template<>
    struct BiquadCascade<q31_t> {
        typedef arm_biquad_casd_df1_inst_q31 type;

        static void init(type *cascade, uint8_t num_stages, q31_t *coefficients, q31_t *delay_pipeline) {
            arm_biquad_cascade_df1_init_q31(cascade, num_stages, coefficients, delay_pipeline, 1);
        }

        static void process_cascade(type *cascade, q31_t *src, q31_t *dst, uint32_t block_size) {
            arm_biquad_cascade_df1_q31(cascade, src, dst, block_size);
        }

        [[nodiscard]] static q31_t multiply(q31_t x, q31_t y) {
            T product;
            arm_mult_q31(&x, &y, &product, 1);
            return product;
        }

        static void
        push_biquad_coefficients(q31_t *coefficients, double b0, double b1, double b2, double a1, double a2) {
            const double max_val = std::max({std::abs(b0), std::abs(b1), std::abs(b2),
                                             std::abs(a1), std::abs(a2)});
            if (max_val > 1.0) {
                // Normalize coefficients if necessary
                double inv_max_val = 1.0 / max_val;
                b0 *= inv_max_val;
                b1 *= inv_max_val;
                b2 *= inv_max_val;
                a1 *= inv_max_val;
                a2 *= inv_max_val;
            }
            const double coefficients_double[] = {b0, b1, b2, -a1, -a2};
            arm_f64_to_q31(coefficients_double, coefficients, COEFFICIENTS_PER_BIQUAD_BLOCK);
        }
    };

    [[nodiscard]] const T *get_coefficients() const {
        return _coefficients;
    }

#if CASCADE_FILTER_DEBUG > 0

    void print_coefficients() const {
        double gain_double;
        to_double(&_gain, &gain_double, 1);
        char debug_buf[128];
        snprintf(debug_buf, sizeof(debug_buf), "gain: %.8f\n", gain_double);
        printf("%s", debug_buf);
        snprintf(debug_buf, sizeof(debug_buf), "coefficients:\n");
        printf("%s", debug_buf);
        for (int i = 0; i < NUMBER_OF_BIQUAD_BLOCKS; ++i) {
            double biquad_coefficients_double[COEFFICIENTS_PER_BIQUAD_BLOCK];
            to_double(&_coefficients[COEFFICIENTS_PER_BIQUAD_BLOCK * i],
                      biquad_coefficients_double, COEFFICIENTS_PER_BIQUAD_BLOCK);
            snprintf(debug_buf, sizeof(debug_buf),
                     "%.8f, %.8f, %.8f, 1, %.8f, %.8f\n",
                     biquad_coefficients_double[0],
                     biquad_coefficients_double[1],
                     biquad_coefficients_double[2],
                     biquad_coefficients_double[3],
                     biquad_coefficients_double[4]);
            printf("%s", debug_buf);
        }
    }

#endif

    bool push_biquad_coefficients(double b0, double b1, double b2, double a1, double a2) {
        if (_num_biquad_blocks_set >= NUMBER_OF_BIQUAD_BLOCKS) {
            return false;
        }
        T *current_coefficients_block = &_coefficients[_num_biquad_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK];
        BiquadCascade<T>::push_biquad_coefficients(current_coefficients_block, b0, b1, b2, a1, a2);
        ++_num_biquad_blocks_set;
        return true;
    }

    void init_biquad_cascades() {
        BiquadCascade<T>::init(&_biquad_cascade, NUMBER_OF_BIQUAD_BLOCKS, _coefficients, _delay_pipeline);
    }

    void reset_state() {
        memset(_delay_pipeline, 0, sizeof(_delay_pipeline));
    }

    void reset() {
        _num_biquad_blocks_set = 0;
        memset(_coefficients, 0, sizeof(_coefficients));
        memset(_delay_pipeline, 0, sizeof(_delay_pipeline));
    }

    [[nodiscard]] T get_gain() const {
        return _gain;
    }

    [[nodiscard]] int get_number_of_blocks() const {
        return _num_biquad_blocks_set;
    }

    void set_gain(double gain) {
        if constexpr (std::is_same<T, double>::value) {
            _gain = gain;
        } else {
            to_native(&gain, &_gain, 1);
        }
    }

    [[nodiscard]] T process(T x) {
        x = BiquadCascade<T>::multiply(x, _gain);
        T output = 0;
        BiquadCascade<T>::process_cascade(&_biquad_cascade, &x, &output, 1);
        return output;
    }

    void process(const T *x, T *out, uint32_t num_samples) {
        if (num_samples == 0) {
            return;
        }
        T input_buffer[num_samples];
        for (uint32_t i = 0; i < num_samples; ++i) {
            input_buffer[i] = BiquadCascade<T>::multiply(x[i], _gain);
        }
        BiquadCascade<T>::process_cascade(&_biquad_cascade, input_buffer, out, num_samples);
    }

    [[nodiscard]] T process(const T *x, uint32_t num_samples) {
        if (num_samples == 0) {
            return 0;
        }
        T output_buffer[num_samples];
        process(x, output_buffer, num_samples);
        return output_buffer[num_samples - 1];
    }

    template<typename U, typename std::enable_if<!std::is_same<T, U>::value, int>::type = 0>
    [[nodiscard]] T process(U x) {
        T x_native;
        to_native(&x, &x_native, 1);
        x_native = BiquadCascade<T>::multiply(x_native, _gain);
        T output = 0;
        BiquadCascade<T>::process_cascade(&_biquad_cascade, &x_native, &output, 1);
        return output;
    }

    template<typename U, typename std::enable_if<!std::is_same<T, U>::value, int>::type = 0>
    void process(const U *x, T *out, uint32_t num_samples) {
        if (num_samples == 0) {
            return;
        }
        T x_native[num_samples];
        to_native(x, x_native, num_samples);
        for (uint32_t i = 0; i < num_samples; ++i) {
            x_native[i] = BiquadCascade<T>::multiply(x_native[i], _gain);
        }
        BiquadCascade<T>::process_cascade(&_biquad_cascade, x_native, out, num_samples);
    }

    template<typename U, typename std::enable_if<!std::is_same<T, U>::value, int>::type = 0>
    [[nodiscard]] T process(const U *x, uint32_t num_samples) {
        if (num_samples == 0) {
            return 0;
        }
        T output_buffer[num_samples];
        process(x, output_buffer, num_samples);
        return output_buffer[num_samples - 1];
    }

    static constexpr unsigned int NUMBER_OF_BIQUAD_BLOCKS = (N + 1) / 2;
    static constexpr unsigned int COEFFICIENTS_PER_BIQUAD_BLOCK = 5;
    static constexpr unsigned int NUMBER_OF_COEFFICIENTS
            = NUMBER_OF_BIQUAD_BLOCKS * COEFFICIENTS_PER_BIQUAD_BLOCK;
    static constexpr unsigned int DELAY_LINE_SIZE = std::is_same<T, q31_t>::value
                                                    ? NUMBER_OF_BIQUAD_BLOCKS * 4
                                                    : NUMBER_OF_BIQUAD_BLOCKS * 2;

private:
    typedef typename BiquadCascade<T>::type biquad_cascade_t;

    biquad_cascade_t _biquad_cascade;

    T _gain{};
    T _coefficients[NUMBER_OF_COEFFICIENTS]{};
    uint32_t _num_biquad_blocks_set = 0;

    T _delay_pipeline[DELAY_LINE_SIZE]{};
};

}
