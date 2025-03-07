#pragma once

#include <dsp/filtering_functions.h>
#include <type_traits>

namespace tiny_iir {

static constexpr unsigned int COEFFICIENTS_PER_BIQUAD_BLOCK = 5;

template<typename>
struct BiquadCascade;

template<>
struct BiquadCascade<float> {
    typedef arm_biquad_cascade_df2T_instance_f32 type;

    static void init(type *cascade,
                     uint8_t num_stages,
                     float *coefficients,
                     float *delay_pipeline) {
        arm_biquad_cascade_df2T_init_f32(
                cascade, num_stages, coefficients, delay_pipeline);
    }

    static void
    process_cascade(type *cascade, float *src, float *dst, uint32_t block_size) {
        arm_biquad_cascade_df2T_f32(cascade, src, dst, block_size);
    }

    [[nodiscard]] static float multiply(float x, float y) {
        return x * y;
    }

    static void push_biquad_coefficients(float *coefficients,
                                         double b0,
                                         double b1,
                                         double b2,
                                         double a1,
                                         double a2) {
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

    static void init(type *cascade,
                     uint8_t num_stages,
                     double *coefficients,
                     double *delay_pipeline) {
        arm_biquad_cascade_df2T_init_f64(
                cascade, num_stages, coefficients, delay_pipeline);
    }

    static void process_cascade(type *cascade,
                                double *src,
                                double *dst,
                                uint32_t block_size) {
        arm_biquad_cascade_df2T_f64(cascade, src, dst, block_size);
    }

    [[nodiscard]] static double multiply(double x, double y) {
        return x * y;
    }

    static void push_biquad_coefficients(double *coefficients,
                                         double b0,
                                         double b1,
                                         double b2,
                                         double a1,
                                         double a2) {
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

    static void init(type *cascade,
                     uint8_t num_stages,
                     q31_t *coefficients,
                     q31_t *delay_pipeline) {
        arm_biquad_cascade_df1_init_q31(
                cascade, num_stages, coefficients, delay_pipeline, 1);
    }

    static void
    process_cascade(type *cascade, q31_t *src, q31_t *dst, uint32_t block_size) {
        arm_biquad_cascade_df1_q31(cascade, src, dst, block_size);
    }

    [[nodiscard]] static q31_t multiply(q31_t x, q31_t y) {
        q31_t product;
        arm_mult_q31(&x, &y, &product, 1);
        return product;
    }

    static void push_biquad_coefficients(q31_t *coefficients,
                                         double b0,
                                         double b1,
                                         double b2,
                                         double a1,
                                         double a2) {
        const double max_val = std::max({std::abs(b0),
                                         std::abs(b1),
                                         std::abs(b2),
                                         std::abs(a1),
                                         std::abs(a2)});
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
        arm_f64_to_q31(
                coefficients_double, coefficients, COEFFICIENTS_PER_BIQUAD_BLOCK);
    }
};

} // namespace tiny_iir