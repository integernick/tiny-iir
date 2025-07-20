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
        // Decide postShift based on the largest absolute coefficient seen across all sections.
        // If no sections were staged yet, fall back to postShift = 0.
        int k = 0;
        if (_max_abs > 1.0) {
            k = static_cast<int>(std::ceil(std::log2(_max_abs)));
            if (k < 0) {
                k = 0;
            }
        }
        _post_shift = static_cast<int8_t>(k);

        // Scale factor 2^{-k}
        const double scale = std::ldexp(1.0, -k);
        const double one_minus = std::nextafter(1.0, 0.0);

        // Write scaled coefficients per section into the caller's buffer.
        for (int s = 0; s < _staged; ++s) {
            const int base = 5 * s;
            double scaled[5];
            for (int i = 0; i < 5; ++i) {
                double x = _raw[base + i] * scale;
                if (x > one_minus) x = one_minus;
                if (x < -1.0) x = -1.0;
                scaled[i] = x;
            }
            arm_f64_to_q31(scaled, &coefficients[base], 5);
        }

        // Initialize CMSIS instance with computed postShift.
        arm_biquad_cascade_df1_init_q31(
                &_cascade_instance,
                num_stages,
                coefficients,
                delay_pipeline,
                _post_shift);

        // reset staging for potential next design
        _staged = 0;
        _max_abs = 0.0;
        _post_shift = 0;
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
    void push_biquad_coefficients(q31_t *, BiquadCoefficients &biquad_coefficients) {
        if (_staged >= MAX_SECTIONS) {
            // Ignore extra sections (or assert in debug builds)
            return;
        }

        // CMSIS-DSP expects a1,a2 with inverted sign in memory.
        const double b0 = biquad_coefficients.b0;
        const double b1 = biquad_coefficients.b1;
        const double b2 = biquad_coefficients.b2;
        const double a1 = -biquad_coefficients.a1;
        const double a2 = -biquad_coefficients.a2;

        const int base = 5 * _staged;
        _raw[base + 0] = b0;
        _raw[base + 1] = b1;
        _raw[base + 2] = b2;
        _raw[base + 3] = a1;
        _raw[base + 4] = a2;

        const double section_max = std::max({std::abs(b0), std::abs(b1), std::abs(b2),
                                             std::abs(a1), std::abs(a2)});
        if (section_max > _max_abs) _max_abs = section_max;

        ++_staged;
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 4;

private:
    arm_biquad_casd_df1_inst_q31 _cascade_instance;
    static constexpr int MAX_SECTIONS = 10; // supports filter order up to 20
    int _staged = 0;
    double _raw[5 * MAX_SECTIONS]{};
    double _max_abs = 0.0;
    int8_t _post_shift = 0;
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
        const double max_val = std::max({std::abs(biquad_coefficients.b0),
                                         std::abs(biquad_coefficients.b1),
                                         std::abs(biquad_coefficients.b2)});
        BiquadCoefficients normalized_biquad_coefficients = {
                .b0 = biquad_coefficients.b0,
                .b1 = biquad_coefficients.b1,
                .b2 = biquad_coefficients.b2,
                .a1 = -biquad_coefficients.a1,
                .a2 = -biquad_coefficients.a2
        };

        if (max_val > 1.0) {
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