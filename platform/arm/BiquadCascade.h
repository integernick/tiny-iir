#pragma once

#include <common/common_utils.h>
#include <dsp/filtering_functions.h>

#include <type_traits>

namespace tiny_iir {

namespace {
struct ScientificNotation {
    double scale = 1.0;
    double mag = 0;
    int8_t exp = 0;
};

ScientificNotation get_scientific_notation(double value) {
    ScientificNotation sn{
            .scale = 1.0,
            .mag = value,
            .exp = 0,
    };

    if (value > 1.0) {
        // Scale the gain inverse to the range [0.5, 1.0)
        auto k = static_cast<int8_t>(std::ceil(std::log2(value)));
        if (k < 0) {
            k = 0;
        }
        sn.exp = k;

        // Scale factor 2^{-k}
        sn.scale = std::ldexp(1.0, -k);
        sn.mag = value * sn.scale;
    }

    return sn;
}

}

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
     * @brief   Get the number of biquad blocks set.
     *
     * @return  The number of biquad blocks set.
     */
    [[nodiscard]] uint32_t get_num_biquad_blocks_set() const {
        return _num_biquad_blocks_set;
    }

    /**
     * @brief   Set the coefficients.
     *
     * @param coefficients  The coefficients array.
     */
    void set_coefficients(float *coefficients) {
        _coefficients = coefficients;
    }

    /**
     * @brief   Initialize the biquad cascade filter.
     *
     * @param num_stages  The number of stages.
     * @param delay_pipeline  The delay pipeline array.
     */
    void init(uint8_t num_stages, float *delay_pipeline) {
        arm_biquad_cascade_df2T_init_f32(
                &_cascade_instance, num_stages, _coefficients, delay_pipeline);
    }

    /**
     * @brief   Reset the filter state.
     */
    void reset() {
        _num_biquad_blocks_set = 0;
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

        _num_biquad_blocks_set++;
    }

    void update_gain(double biquad_gain_inv, float &gain) {
        gain *= biquad_gain_inv;
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 2;

    static constexpr double UNITY = 1.0;
private:
    arm_biquad_cascade_df2T_instance_f32 _cascade_instance;
    float *_coefficients = nullptr;
    uint32_t _num_biquad_blocks_set = 0;
};

/**
 * @brief   Biquad cascade for double native type.
 */
template<>
class BiquadCascade<double> {
public:
    /**
     * @brief   Get the number of biquad blocks set.
     *
     * @return  The number of biquad blocks set.
     */
    [[nodiscard]] uint32_t get_num_biquad_blocks_set() const {
        return _num_biquad_blocks_set;
    }

    /**
     * @brief   Set the coefficients.
     *
     * @param coefficients  The coefficients array.
     */
    void set_coefficients(double *coefficients) {
        _coefficients = coefficients;
    }

    /**
     * @brief   Initialize the biquad cascade filter.
     *
     * @param num_stages  The number of stages.
     * @param delay_pipeline  The delay pipeline array.
     */
    void init(uint8_t num_stages, double *delay_pipeline) {
        arm_biquad_cascade_df2T_init_f64(
                &_cascade_instance, num_stages, _coefficients, delay_pipeline);
    }

    /**
     * @brief   Reset the filter state.
     */
    void reset() {
        _num_biquad_blocks_set = 0;
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
        coefficients[0] = biquad_coefficients.b0;
        coefficients[1] = biquad_coefficients.b1;
        coefficients[2] = biquad_coefficients.b2;
        coefficients[3] = -biquad_coefficients.a1;
        coefficients[4] = -biquad_coefficients.a2;

        _num_biquad_blocks_set++;
    }

    void update_gain(double biquad_gain_inv, double &gain) {
        gain *= biquad_gain_inv;
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 2;

    static constexpr float UNITY = 1.0f;
private:
    arm_biquad_cascade_df2T_instance_f64 _cascade_instance;
    double *_coefficients = nullptr;
    uint32_t _num_biquad_blocks_set = 0;
};

/**
 * @brief   Biquad cascade for Q31 native type.
 */
template<>
class BiquadCascade<q31_t> {
public:
    /**
     * @brief   Get the number of biquad blocks set.
     *
     * @return  The number of biquad blocks set.
     */
    [[nodiscard]] uint32_t get_num_biquad_blocks_set() const {
        return _num_biquad_blocks_set;
    }

    /**
     * @brief   Set the coefficients.
     *
     * @param coefficients  The coefficients array.
     */
    void set_coefficients(q31_t *coefficients) {
        _coefficients = coefficients;
    }

    /**
     * @brief   Initialize the biquad cascade filter.
     *
     * @param num_stages  The number of stages.
     * @param delay_pipeline  The delay pipeline array.
     */
    void init(uint8_t num_stages, q31_t *delay_pipeline) {
        // Initialize CMSIS instance with computed postShift.
        arm_biquad_cascade_df1_init_q31(
                &_cascade_instance,
                num_stages,
                _coefficients,
                delay_pipeline,
                _post_shift);
    }

    /**
     * @brief   Reset the filter state.
     */
    void reset() {
        _post_shift = 0;
        _num_biquad_blocks_set = 0;
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
        arm_shift_q31(dst, _gain_shift, dst, block_size);
    }

    /**
     * @brief   Push biquad coefficients to the filter.
     *
     * @param coefficients  The pointer to the biquad coefficients block.
     * @param biquad_coefficients  The biquad coefficients.
     */
    void push_biquad_coefficients(q31_t *coefficients, const BiquadCoefficients &biquad_coefficients) {
        BiquadCoefficients normalized_biquad_coefficients = {
                .b0 = biquad_coefficients.b0,
                .b1 = biquad_coefficients.b1,
                .b2 = biquad_coefficients.b2,
                .a1 = -biquad_coefficients.a1,
                .a2 = -biquad_coefficients.a2
        };

        const double max_val = std::max({std::abs(biquad_coefficients.b0),
                                         std::abs(biquad_coefficients.b1),
                                         std::abs(biquad_coefficients.b2),
                                         std::abs(biquad_coefficients.a1),
                                         std::abs(biquad_coefficients.a2)});

        double scale = 1.0;
        if (max_val > 1.0) {
            // Normalize coefficients if necessary
            const ScientificNotation sn = get_scientific_notation(max_val);
            scale = sn.scale;
            update_post_shift(sn.exp);
        } else {
            // Each biquad block is shifted by _post_shift bits, so the coefficients need to be scaled in any case
            scale = 1.0 / static_cast<double>(1UL << _post_shift);
        }

        normalized_biquad_coefficients.b0 *= scale;
        normalized_biquad_coefficients.b1 *= scale;
        normalized_biquad_coefficients.b2 *= scale;
        normalized_biquad_coefficients.a1 *= scale;
        normalized_biquad_coefficients.a2 *= scale;


        arm_f64_to_q31(&normalized_biquad_coefficients.b0, coefficients, COEFFICIENTS_PER_BIQUAD_BLOCK);

        _num_biquad_blocks_set++;

        double coeffs_double[COEFFICIENTS_PER_BIQUAD_BLOCK];
        to_double(coefficients, coeffs_double, COEFFICIENTS_PER_BIQUAD_BLOCK);
    }

    void update_gain(double biquad_gain_inv, q31_t &gain) {
        /*const ScientificNotation sn = get_scientific_notation(biquad_gain_inv);
        update_post_shift(sn.exp);
        q31_t biquad_gain_inv_q31;
        arm_f64_to_q31(&sn.mag, &biquad_gain_inv_q31, 1);
        multiply(&gain, &biquad_gain_inv_q31, &gain, 1);
*/

        double g = biquad_gain_inv;
        int8_t e = 0;
        if (g > 1.0) {
            e = (int8_t)std::ceil(std::log2(g));
            if (e < 0) e = 0;
            if (e > 30) e = 30;                      // sanity cap
            g = std::ldexp(g, -e);                   // g := g / 2^e  in (0,1]
        }
        // Clamp mantissa away from 1.0 to avoid rounding to 2^31 in Q31
        const double one_minus = std::nextafter(1.0, 0.0);
        if (g > one_minus) g = one_minus;

        q31_t m_q31;
        arm_f64_to_q31(&g, &m_q31, 1);               // m in Q31

        // Multiply current Q31 gain by mantissa
        q31_t prod;
        arm_mult_q31(&gain, &m_q31, &prod, 1);       // still <= 1.0 in Q31

        // Try to absorb as much of 2^e as possible into the mantissa by left-shifting prod
        int8_t hr = headroom_q31(prod);
        int8_t absorb = (e < hr) ? e : hr;
        if (absorb > 0) {
            prod = (q31_t)((uint32_t)prod << absorb);
            e -= absorb;
        }

        // Commit results
        gain = prod;                                  // new mantissa
        _gain_shift = (int8_t)(_gain_shift + e);      // leftover power-of-two (global, not per-stage)
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 4;

    static constexpr q31_t UNITY = Q31_MAX;

private:
    static inline int8_t headroom_q31(q31_t x) {
        if (x == 0) return 30;                       // plenty of room
        uint32_t ax = (x < 0) ? (uint32_t)(-x) : (uint32_t)x;
        // Highest allowed magnitude bit is at position 30 (0x4000_0000).
        // headroom = number of shifts before hitting bit30.
        int lz = __builtin_clz(ax);                  // 0..32
        int hr = lz - 1;                             // keep sign bit at 31, value MSB at <=30
        return (hr < 0) ? 0 : (hr > 30 ? 30 : hr);
    }

    void update_post_shift(int8_t new_shift) {
        if (new_shift > _post_shift) {
            const int8_t delta_shift = new_shift - _post_shift;
            _post_shift = new_shift;
            const uint32_t already_set_block_size = _num_biquad_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK;
            arm_shift_q31(_coefficients, -delta_shift, _coefficients, already_set_block_size);
        }
    }

    arm_biquad_casd_df1_inst_q31 _cascade_instance;
    q31_t *_coefficients = nullptr;
    uint32_t _num_biquad_blocks_set = 0;
    int8_t _post_shift = 0;
    int8_t _gain_shift = 0;
};

/**
 * @brief   Biquad cascade for Q15 native type.
 */
template<>
class BiquadCascade<q15_t> {
public:
    /**
     * @brief   Initialize the biquad cascade filter.
     *
     * @param num_stages  The number of stages.
     * @param delay_pipeline  The delay pipeline array.
     */
    void init(uint8_t num_stages,
              q15_t *delay_pipeline) {
        arm_biquad_cascade_df1_init_q15(
                &_cascade_instance, num_stages, _coefficients, delay_pipeline,
                _post_shift);
    }

    /**
     * @brief   Reset the filter state.
     */
    void reset() {
        _post_shift = 0;
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
        BiquadCoefficients normalized_biquad_coefficients = {
                .b0 = biquad_coefficients.b0,
                .b1 = biquad_coefficients.b1,
                .b2 = biquad_coefficients.b2,
                .a1 = -biquad_coefficients.a1,
                .a2 = -biquad_coefficients.a2
        };

        const double max_val = std::max({std::abs(biquad_coefficients.b0),
                                         std::abs(biquad_coefficients.b1),
                                         std::abs(biquad_coefficients.b2),
                                         std::abs(biquad_coefficients.a1),
                                         std::abs(biquad_coefficients.a2)});

        if (max_val > 1.0) {
            // Normalize coefficients if necessary
            double inv_max_val = 1.0 / max_val;
            normalized_biquad_coefficients.b0 *= inv_max_val;
            normalized_biquad_coefficients.b1 *= inv_max_val;
            normalized_biquad_coefficients.b2 *= inv_max_val;
            normalized_biquad_coefficients.a1 *= inv_max_val;
            normalized_biquad_coefficients.a2 *= inv_max_val;
        }

        arm_f64_to_q15(&normalized_biquad_coefficients.b0, coefficients, COEFFICIENTS_PER_BIQUAD_BLOCK);
    }

    void update_gain(double biquad_gain_inv, q15_t &gain) {
        if (biquad_gain_inv > 1.0) {
            // Scale the gain inverse to the range [0.5, 1.0)
            auto k = static_cast<int8_t>(std::ceil(std::log2(biquad_gain_inv)));
            if (k < 0) {
                k = 0;
            }
            _post_shift += k;

            // Scale factor 2^{-k}
            const double scale = std::ldexp(1.0, -k);
            const double one_minus = std::nextafter(1.0, 0.0);

            // Write scaled coefficients per section into the caller's buffer.
            biquad_gain_inv *= scale;
            if (biquad_gain_inv > one_minus) biquad_gain_inv = one_minus;
            if (biquad_gain_inv < -1.0) biquad_gain_inv = -1.0;
        }

        q15_t biquad_gain_inv_q15;
        arm_f64_to_q15(&biquad_gain_inv, &biquad_gain_inv_q15, 1);
        multiply(&gain, &biquad_gain_inv_q15, &gain, 1);
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 4;

    static constexpr q15_t UNITY = Q15_MAX;
private:
    arm_biquad_casd_df1_inst_q15 _cascade_instance;
    q15_t *_coefficients = nullptr;
    int8_t _post_shift = 0;
    uint32_t _num_coefficients = 0;
};

} // namespace tiny_iir