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
class BiquadCascade<float> {
public:
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
     * @brief   Get the normalizing gain.
     *
     * @return  The normalizing gain.
     */
    [[nodiscard]] double get_normalizing_gain() const {
        return 1.0;
    }

    /**
     * @brief   Initialize the biquad cascade filter.
     *
     * @param num_stages  The number of stages.
     * @param coefficients  The coefficients array.
     * @param delay_pipeline  The delay pipeline array.
     */
    void init(uint8_t num_stages,
              float *coefficients,
              float *delay_pipeline) {
        arm_biquad_cascade_df2T_init_f32(
                &_cascade_instance, num_stages, coefficients, delay_pipeline);
    }

    /**
     * @brief   Process a batch of samples.
     *
     * @param src  The input buffer.
     * @param dst  The output buffer.
     * @param block_size  The number of samples to process.
     */
    void process_cascade(float *src, float *dst, uint32_t block_size) {
        arm_biquad_cascade_df2T_f32(&_cascade_instance, src, dst, block_size);
    }

    /**
     * @brief   Push biquad coefficients to the filter.
     *
     * @param coefficients  The pointer to the biquad coefficients block.
     * @param biquad_coefficients  The biquad coefficients.
     */
    void push_biquad_coefficients(float *coefficients, const BiquadCoefficients &biquad_coefficients) {
        coefficients[0] = static_cast<float>(biquad_coefficients.b0);
        coefficients[1] = static_cast<float>(biquad_coefficients.b1);
        coefficients[2] = static_cast<float>(biquad_coefficients.b2);
        coefficients[3] = static_cast<float>(-biquad_coefficients.a1);
        coefficients[4] = static_cast<float>(-biquad_coefficients.a2);
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 2;

private:
    arm_biquad_cascade_df2T_instance_f32 _cascade_instance;
};

/**
 * @brief   Biquad cascade for double native type.
 */
template<>
class BiquadCascade<double> {
public:
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
     * @brief   Get the normalizing gain.
     *
     * @return  The normalizing gain.
     */
    [[nodiscard]] double get_normalizing_gain() const {
        return _normalizing_gain;
    }

    /**
     * @brief   Initialize the biquad cascade filter.
     *
     * @param num_stages  The number of stages.
     * @param coefficients  The coefficients array.
     * @param delay_pipeline  The delay pipeline array.
     */
    void init(uint8_t num_stages, double *coefficients, double *delay_pipeline) {
        arm_biquad_cascade_df2T_init_f64(
                &_cascade_instance, num_stages, coefficients, delay_pipeline);
    }

    /**
     * @brief   Process a batch of samples.
     *
     * @param src  The input buffer.
     * @param dst  The output buffer.
     * @param block_size  The number of samples to process.
     */
    void process_cascade(double *src, double *dst, uint32_t block_size) {
        arm_biquad_cascade_df2T_f64(&_cascade_instance, src, dst, block_size);
    }

    /**
     * @brief   Push biquad coefficients to the filter.
     *
     * @param coefficients  The pointer to the biquad coefficients block.
     * @param biquad_coefficients  The biquad coefficients.
     */
    void push_biquad_coefficients(double *coefficients, const BiquadCoefficients &biquad_coefficients) {
         double max_val = std::max({std::abs(biquad_coefficients.b0),
                                         std::abs(biquad_coefficients.b1),
                                         std::abs(biquad_coefficients.b2)});
        BiquadCoefficients normalized_biquad_coefficients = {
                .b0 = biquad_coefficients.b0,
                .b1 = biquad_coefficients.b1,
                .b2 = biquad_coefficients.b2,
                .a1 = -biquad_coefficients.a1,
                .a2 = -biquad_coefficients.a2
        };

        coefficients[0] = normalized_biquad_coefficients.b0;
        coefficients[1] = normalized_biquad_coefficients.b1;
        coefficients[2] = normalized_biquad_coefficients.b2;
        coefficients[3] = normalized_biquad_coefficients.a1;
        coefficients[4] = normalized_biquad_coefficients.a2;
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 2;

private:
    double _normalizing_gain = 1.0;
    arm_biquad_cascade_df2T_instance_f64 _cascade_instance;
};

/**
 * @brief   Biquad cascade for Q31 native type.
 */
template<>
class BiquadCascade<q31_t> {
public:
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

    /**sd
     * @brief   Initialize the biquad cascade filter.
     *
     * @param num_stages  The number of stages.
     * @param coefficients  The coefficients array.
     * @param delay_pipeline  The delay pipeline array.
     */
    void init(uint8_t num_stages,
              q31_t *coefficients,
              q31_t *delay_pipeline) {
        arm_biquad_cascade_df1_init_q31(
                &_cascade_instance, num_stages, coefficients, delay_pipeline,
                0);
    }

    /**
     * @brief   Process a batch of samples.
     *
     * @param src  The input buffer.
     * @param dst  The output buffer.
     * @param block_size  The number of samples to process.
     */
    void process_cascade(q31_t *src, q31_t *dst, uint32_t block_size) {
        arm_biquad_cascade_df1_q31(&_cascade_instance, src, dst, block_size);
    }

    /**
     * @brief   Push biquad coefficients to the filter.
     *
     * @param coefficients  The pointer to the biquad coefficients block.
     * @param biquad_coefficients  The biquad coefficients.
     */
    void push_biquad_coefficients(q31_t *coefficients, BiquadCoefficients &biquad_coefficients) {
        double max_val = std::max({std::abs(biquad_coefficients.b0),
                                         std::abs(biquad_coefficients.b1),
                                         std::abs(biquad_coefficients.b2)});

        if (max_val > 1.0) {
            // Normalize coefficients if necessary
            const double scale = 1.0 / max_val;

            biquad_coefficients.b0 *= scale;
            biquad_coefficients.b1 *= scale;
            biquad_coefficients.b2 *= scale;
        }

        /* CMSIS-DSP a1, a2 signs are inverted */
        biquad_coefficients.a1 = -biquad_coefficients.a1;
        biquad_coefficients.a2 = -biquad_coefficients.a2;
        arm_f64_to_q31(&biquad_coefficients.b0, coefficients, COEFFICIENTS_PER_BIQUAD_BLOCK);
        biquad_coefficients.a1 = -biquad_coefficients.a1;
        biquad_coefficients.a2 = -biquad_coefficients.a2;
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 4;

private:
    arm_biquad_casd_df1_inst_q31 _cascade_instance;
};

/**
 * @brief   Biquad cascade for Q31 native type.
 */
template<>
class BiquadCascade<q15_t> {
public:
    /**
     * @brief   Multiply two numbers (Q31 version).
     *
     * @param x  The first number.
     * @param y  The second number.
     * @return  The product of the two numbers.
     */
    [[nodiscard]] static q15_t multiply(q15_t x, q15_t y) {
        q15_t product;
        arm_mult_q15(&x, &y, &product, 1);
        return product;
    }

    /**sd
     * @brief   Initialize the biquad cascade filter.
     *
     * @param num_stages  The number of stages.
     * @param coefficients  The coefficients array.
     * @param delay_pipeline  The delay pipeline array.
     */
    void init(uint8_t num_stages,
              q15_t *coefficients,
              q15_t *delay_pipeline) {
        arm_biquad_cascade_df1_init_q15(
                &_cascade_instance, num_stages, coefficients, delay_pipeline,
                0);
    }

    /**
     * @brief   Process a batch of samples.
     *
     * @param src  The input buffer.
     * @param dst  The output buffer.
     * @param block_size  The number of samples to process.
     */
    void process_cascade(q15_t *src, q15_t *dst, uint32_t block_size) {
        arm_biquad_cascade_df1_q15(&_cascade_instance, src, dst, block_size);
    }

    /**
     * @brief   Push biquad coefficients to the filter.
     *
     * @param coefficients  The pointer to the biquad coefficients block.
     * @param biquad_coefficients  The biquad coefficients.
     */
    void push_biquad_coefficients(q15_t *coefficients, const BiquadCoefficients &biquad_coefficients) {
        double max_val = std::max({std::abs(biquad_coefficients.b0),
                                   std::abs(biquad_coefficients.b1),
                                   std::abs(biquad_coefficients.b2)});
        BiquadCoefficients normalized_biquad_coefficients = {
                .b0 = biquad_coefficients.b0,
                .b1 = biquad_coefficients.b1,
                .b2 = biquad_coefficients.b2,
                .a1 = -biquad_coefficients.a1,
                .a2 = -biquad_coefficients.a2
        };

        if (false && max_val > 1.0) {
            // Normalize coefficients if necessary
            const double scale = 1.0 / max_val;

            normalized_biquad_coefficients.b0 *= scale;
            normalized_biquad_coefficients.b1 *= scale;
            normalized_biquad_coefficients.b2 *= scale;
        }

        arm_f64_to_q15(&normalized_biquad_coefficients.b0, coefficients, COEFFICIENTS_PER_BIQUAD_BLOCK);
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 4;

private:
    arm_biquad_casd_df1_inst_q15 _cascade_instance;
};

} // namespace tiny_iir