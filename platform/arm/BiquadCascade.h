#pragma once

#include <common/common_utils.h>
#include <dsp/filtering_functions.h>

#include <type_traits>

namespace tiny_iir {

static constexpr unsigned int COEFFICIENTS_PER_BIQUAD_BLOCK = 5;

template<typename>
struct BiquadCascade;

/**
 * @brief   Biquad cascade for float native type.
 */
template<>
struct BiquadCascade<float> {
    typedef arm_biquad_cascade_df2T_instance_f32 type;

    /**
     * @brief   Initialize the biquad cascade filter.
     *
     * @param cascade  The filter instance.
     * @param num_stages  The number of stages.
     * @param coefficients  The coefficients array.
     * @param delay_pipeline  The delay pipeline array.
     */
    static void init(type *cascade,
                     uint8_t num_stages,
                     float *coefficients,
                     float *delay_pipeline) {
        arm_biquad_cascade_df2T_init_f32(
                cascade, num_stages, coefficients, delay_pipeline);
    }

    /**
     * @brief   Process a batch of samples.
     *
     * @param cascade  The filter instance.
     * @param src  The input buffer.
     * @param dst  The output buffer.
     * @param block_size  The number of samples to process.
     */
    static void process_cascade(type *cascade, float *src, float *dst, uint32_t block_size) {
        arm_biquad_cascade_df2T_f32(cascade, src, dst, block_size);
    }

    /**
     * @brief   Multiply two numbers (float version).
     *
     * @param x  The first number.
     * @param y  The second number.
     * @return  The product of the two numbers.
     */
    [[nodiscard]] static float multiply(float x, float y) {
        return x * y;
    }

    /**
     * @brief   Push biquad coefficients to the filter.
     *
     * @param coefficients  The pointer to the biquad coefficients block.
     * @param biquad_coefficients  The biquad coefficients.
     */
    static void push_biquad_coefficients(float *coefficients, const BiquadCoefficients &biquad_coefficients) {
        coefficients[0] = static_cast<float>(biquad_coefficients.b0);
        coefficients[1] = static_cast<float>(biquad_coefficients.b1);
        coefficients[2] = static_cast<float>(biquad_coefficients.b2);
        coefficients[3] = static_cast<float>(-biquad_coefficients.a1);
        coefficients[4] = static_cast<float>(-biquad_coefficients.a2);
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 2;
};

/**
 * @brief   Biquad cascade for double native type.
 */
template<>
struct BiquadCascade<double> {
    typedef arm_biquad_cascade_df2T_instance_f64 type;

    /**
     * @brief   Initialize the biquad cascade filter.
     *
     * @param cascade  The filter instance.
     * @param num_stages  The number of stages.
     * @param coefficients  The coefficients array.
     * @param delay_pipeline  The delay pipeline array.
     */
    static void init(type *cascade,
                     uint8_t num_stages,
                     double *coefficients,
                     double *delay_pipeline) {
        arm_biquad_cascade_df2T_init_f64(
                cascade, num_stages, coefficients, delay_pipeline);
    }

    /**
     * @brief   Process a batch of samples.
     *
     * @param cascade  The filter instance.
     * @param src  The input buffer.
     * @param dst  The output buffer.
     * @param block_size  The number of samples to process.
     */
    static void process_cascade(type *cascade,
                                double *src,
                                double *dst,
                                uint32_t block_size) {
        arm_biquad_cascade_df2T_f64(cascade, src, dst, block_size);
    }

    /**
     * @brief   Multiply two numbers (double version).
     *
     * @param x  The first number.
     * @param y  The second number.
     * @return  The product of the two numbers.
     */
    [[nodiscard]] static double multiply(double x, double y) {
        return x * y;
    }

    /**
     * @brief   Push biquad coefficients to the filter.
     *
     * @param coefficients  The pointer to the biquad coefficients block.
     * @param biquad_coefficients  The biquad coefficients.
     */
    static void push_biquad_coefficients(double *coefficients, const BiquadCoefficients &biquad_coefficients) {
        coefficients[0] = biquad_coefficients.b0;
        coefficients[1] = biquad_coefficients.b1;
        coefficients[2] = biquad_coefficients.b2;
        coefficients[3] = -biquad_coefficients.a1;
        coefficients[4] = -biquad_coefficients.a2;
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 2;
};

/**
 * @brief   Biquad cascade for Q31 native type.
 */
template<>
struct BiquadCascade<q31_t> {

    typedef arm_biquad_casd_df1_inst_q31 type;

    /**
     * @brief   Initialize the biquad cascade filter.
     *
     * @param cascade  The filter instance.
     * @param num_stages  The number of stages.
     * @param coefficients  The coefficients array.
     * @param delay_pipeline  The delay pipeline array.
     */
    static void init(type *cascade,
                     uint8_t num_stages,
                     q31_t *coefficients,
                     q31_t *delay_pipeline) {
        arm_biquad_cascade_df1_init_q31(
                cascade, num_stages, coefficients, delay_pipeline, 1);
    }

    /**
     * @brief   Process a batch of samples.
     *
     * @param cascade  The filter instance.
     * @param src  The input buffer.
     * @param dst  The output buffer.
     * @param block_size  The number of samples to process.
     */
    static void process_cascade(type *cascade, q31_t *src, q31_t *dst, uint32_t block_size) {
        arm_biquad_cascade_df1_q31(cascade, src, dst, block_size);
    }

    /**
     * @brief   Multiply two numbers (Q31 version).
     *
     * @param x  The first number.
     * @param y  The second number.
     * @return  The product of the two numbers.
     */
    [[nodiscard]] static q31_t multiply(q31_t x, q31_t y) {
        q31_t product;
        arm_mult_q31(&x, &y, &product, 1);
        return product;
    }

    /**
     * @brief   Push biquad coefficients to the filter.
     *
     * @param coefficients  The pointer to the biquad coefficients block.
     * @param biquad_coefficients  The biquad coefficients.
     */
    static void push_biquad_coefficients(q31_t *coefficients, const BiquadCoefficients &biquad_coefficients) {
        const double max_val = std::max({std::abs(biquad_coefficients.b0),
                                         std::abs(biquad_coefficients.b1),
                                         std::abs(biquad_coefficients.b2),
                                         std::abs(biquad_coefficients.a1),
                                         std::abs(biquad_coefficients.a2)});
        BiquadCoefficients normalized_biquad_coefficients = {
                .b0 = biquad_coefficients.b0,
                .b1 = biquad_coefficients.b1,
                .b2 = biquad_coefficients.b2,
                .a1 = -biquad_coefficients.a1,
                .a2 = -biquad_coefficients.a2
        };

        if (max_val > 1.0) {
            // Normalize coefficients if necessary
            double inv_max_val = 1.0 / max_val;
            normalized_biquad_coefficients.b0 *= inv_max_val;
            normalized_biquad_coefficients.b1 *= inv_max_val;
            normalized_biquad_coefficients.b2 *= inv_max_val;
            normalized_biquad_coefficients.a1 *= inv_max_val;
            normalized_biquad_coefficients.a2 *= inv_max_val;
        }

        arm_f64_to_q31(&normalized_biquad_coefficients.b0, coefficients, COEFFICIENTS_PER_BIQUAD_BLOCK);
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 4;
};

} // namespace tiny_iir