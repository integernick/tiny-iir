#pragma once

#include <common/common_utils.h>
#include <type_utils.h>

#include <dsp/filtering_functions.h>

#include <type_traits>

namespace tiny_iir {

namespace {
struct ScientificNotation {
    double scale = 1.0;
    int8_t exp = 0;
};

ScientificNotation get_scientific_notation(double value) {
    ScientificNotation sn{
            .scale = 1.0,
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
    }

    return sn;
}
} // namespace

template<typename T>
struct coeffs_per_stage {
    // arm_biquad_casd_df1_inst_q31 expects [b0, b1, b2, a1, a2]
    static constexpr uint32_t value = 5;
};
template<>
struct coeffs_per_stage<q15_t> {
    // arm_biquad_casd_df1_inst_q15 expects [b0, 0, b1, b2, a1, a2]
    static constexpr uint32_t value = 6;
};

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

    /**
     * @brief   Get the biquad coefficients (double representation).
     *
     * @param biquad_idx  The index of the biquad.
     * @return  The biquad coefficients.
     */
    [[nodiscard]] BiquadCoefficients get_biquad_coefficients(uint32_t biquad_idx) const {
        BiquadCoefficients coeffs;
        to_double(&_coefficients[biquad_idx * coeffs_per_stage<float>::value], &coeffs.b0,
                  coeffs_per_stage<float>::value);
        coeffs.a1 = -coeffs.a1;
        coeffs.a2 = -coeffs.a2;

        return coeffs;
    }

    /**
     * @brief   Update the cascade gain from the biquad gain inverse.
     *
     * @param biquad_gain_inv  The biquad gain inverse.
     * @param gain  The cascade gain.
     */
    static void update_gain(double biquad_gain_inv, float &gain) {
        gain *= static_cast<float>(biquad_gain_inv);
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 2;

    static constexpr double UNITY = 1.0;

private:
    arm_biquad_cascade_df2T_instance_f32 _cascade_instance{};
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

    /**
     * @brief   Get the biquad coefficients (double representation).
     *
     * @param biquad_idx  The index of the biquad.
     * @return  The biquad coefficients.
     */
    [[nodiscard]] BiquadCoefficients get_biquad_coefficients(uint32_t biquad_idx) const {
        BiquadCoefficients coeffs;
        to_double(&_coefficients[biquad_idx * coeffs_per_stage<double>::value], &coeffs.b0,
                  coeffs_per_stage<double>::value);
        coeffs.a1 = -coeffs.a1;
        coeffs.a2 = -coeffs.a2;

        return coeffs;
    }

    /**
     * @brief   Update the cascade gain from the biquad gain inverse.
     *
     * @param biquad_gain_inv  The biquad gain inverse.
     * @param gain  The cascade gain.
     */
    static void update_gain(double biquad_gain_inv, double &gain) {
        gain *= biquad_gain_inv;
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 2;
    static constexpr float UNITY = 1.0f;

private:
    arm_biquad_cascade_df2T_instance_f64 _cascade_instance{};
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
        _gain_shift = 0;
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

        double scale;

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


        arm_f64_to_q31(&normalized_biquad_coefficients.b0, coefficients, coeffs_per_stage<q31_t>::value);

        _num_biquad_blocks_set++;
    }

    /**
     * @brief   Get the biquad coefficients (double representation).
     *
     * @param biquad_idx  The index of the biquad.
     * @return  The biquad coefficients.
     */
    [[nodiscard]] BiquadCoefficients get_biquad_coefficients(uint32_t biquad_idx) const {
        BiquadCoefficients coeffs;
        to_double(&_coefficients[biquad_idx * coeffs_per_stage<q31_t>::value], &coeffs.b0,
                  coeffs_per_stage<q31_t>::value);
        coeffs.a1 = -coeffs.a1;
        coeffs.a2 = -coeffs.a2;

        return coeffs;
    }

    /**
     * @brief   Update the cascade gain from the biquad gain inverse.
     *
     * @param biquad_gain_inv  The biquad gain inverse.
     * @param gain  The cascade gain.
     */
    void update_gain(double biquad_gain_inv, q31_t &gain) {
        double g = biquad_gain_inv;
        int8_t e = 0;

        if (g > 1.0) {
            e = (int8_t) std::ceil(std::log2(g));
            e = constrain(e, static_cast<int8_t>(0), static_cast<int8_t>(30));
            g = std::ldexp(g, -e); // g := g / 2^e  in (0,1]
        }

        // Clamp mantissa away from 1.0 to avoid rounding to 2^31 in Q31
        const double one_minus = std::nextafter(1.0, 0.0);
        if (g > one_minus) {
            g = one_minus;
        }

        q31_t m_q31;
        arm_f64_to_q31(&g, &m_q31, 1);

        // Multiply current Q31 gain by mantissa
        q31_t prod;
        arm_mult_q31(&gain, &m_q31, &prod, 1);

        // Try to absorb as much of 2^e as possible into the mantissa by left-shifting prod
        int8_t hr = headroom_q31(prod);
        int8_t absorb = (e < hr) ? e : hr;
        if (absorb > 0) {
            prod = (q31_t) ((uint32_t) prod << absorb);
            e -= absorb;
        }

        gain = prod;
        _gain_shift = (int8_t) (_gain_shift + e);
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 4;
    static constexpr q31_t UNITY = Q31_MAX;

private:
    /**
     * @brief   Computes the number of bits that it is safe to shift the input value left.
     *
     * @param x The q15_t input value.
     * @return  The headroom of the input value.
     */
    static inline int8_t headroom_q31(q31_t x) {
        if (x == 0) {
            return 30;
        }

        const uint32_t ax = (x < 0) ? static_cast<uint32_t>(-x) : static_cast<uint32_t>(x);
        int lz = __builtin_clz(ax);
        int hr = lz - 1;
        hr = constrain(hr, 0, 30);
        return static_cast<int8_t>(hr);
    }

    /**
     * @brief   Update the saved post-shift value.
     *
     * @param new_shift The new post-shift value.
     */
    void update_post_shift(int8_t new_shift) {
        if (new_shift > _post_shift) {
            const int8_t delta_shift = new_shift - _post_shift;
            _post_shift = new_shift;
            const uint32_t already_set_block_size = _num_biquad_blocks_set * coeffs_per_stage<q31_t>::value;
            arm_shift_q31(_coefficients, static_cast<int8_t>(-delta_shift),
                          _coefficients, already_set_block_size);
        }
    }

    arm_biquad_casd_df1_inst_q31 _cascade_instance{};
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
    void set_coefficients(q15_t *coefficients) {
        _coefficients = coefficients;
    }

    /**
     * @brief   Initialize the biquad cascade filter.
     *
     * @param num_stages  The number of stages.
     * @param delay_pipeline  The delay pipeline array.
     */
    void init(uint8_t num_stages, q15_t *delay_pipeline) {
        // Initialize CMSIS instance with computed postShift.
        arm_biquad_cascade_df1_init_q15(
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
        _gain_shift = 0;
        _num_biquad_blocks_set = 0;
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
        arm_shift_q15(dst, _gain_shift, dst, block_size);
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
        double scale;

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

        // arm_biquad_casd_df1_inst_q15 expects [b0, 0, b1, b2, a1, a2]
        arm_f64_to_q15(&normalized_biquad_coefficients.b0, coefficients, 1); // b0
        coefficients[1] = 0;
        arm_f64_to_q15(&normalized_biquad_coefficients.b1, coefficients + 2, 4); // b1, b2, a1, a2

        _num_biquad_blocks_set++;
    }

    /**
     * @brief   Get the biquad coefficients (double representation).
     *
     * @param biquad_idx  The index of the biquad.
     * @return  The biquad coefficients.
     */
    [[nodiscard]] BiquadCoefficients get_biquad_coefficients(uint32_t biquad_idx) const {
        BiquadCoefficients coeffs;
        // arm_biquad_casd_df1_inst_q15 holds [b0, 0, b1, b2, a1, a2]
        to_double(&_coefficients[biquad_idx * coeffs_per_stage<q15_t>::value], &coeffs.b0, 1);
        to_double(&_coefficients[biquad_idx * coeffs_per_stage<q15_t>::value + 2], &coeffs.b1, 4);
        coeffs.a1 = -coeffs.a1;
        coeffs.a2 = -coeffs.a2;

        return coeffs;
    }

    /**
     * @brief   Update the cascade gain from the biquad gain inverse.
     *
     * @param biquad_gain_inv  The biquad gain inverse.
     * @param gain  The cascade gain.
     */
    void update_gain(double biquad_gain_inv, q15_t &gain) {
        double g = biquad_gain_inv;
        int8_t e = 0;

        if (g > 1.0) {
            e = static_cast<int8_t>(std::ceil(std::log2(g)));
            e = constrain(e, static_cast<int8_t>(0), static_cast<int8_t>(14));
            g = std::ldexp(g, -e); // g := g / 2^e  in (0,1]
        }

        // Clamp mantissa away from 1.0 to avoid rounding to 2^15 in Q15
        const double one_minus = std::nextafter(1.0, 0.0);
        if (g > one_minus) {
            g = one_minus;
        }

        q15_t m_q15;
        arm_f64_to_q15(&g, &m_q15, 1);

        // Multiply current Q15 gain by mantissa
        q15_t prod;
        arm_mult_q15(&gain, &m_q15, &prod, 1);

        // Try to absorb as much of 2^e as possible into the mantissa by left-shifting prod
        int8_t hr = headroom_q15(prod);
        int8_t absorb = (e < hr) ? e : hr;
        if (absorb > 0) {
            prod = static_cast<q15_t>((uint32_t) prod << absorb);
            e -= absorb;
        }

        gain = prod;
        _gain_shift = static_cast<int8_t>(_gain_shift + e);
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 4;
    static constexpr q15_t UNITY = Q15_MAX;

private:
    /**
     * @brief   Computes the number of bits that it is safe to shift the input value left.
     *
     * @param x The q15_t input value.
     * @return  The headroom of the input value.
     */
    static inline int8_t headroom_q15(q15_t x) {
        if (x == 0) {
            return 14; // 1 sign + 15 frac -> max safe left shifts is 14
        }

        const uint32_t ax = (x < 0)
                            ? static_cast<uint32_t>(-static_cast<int32_t>(x))
                            : static_cast<uint32_t>(x);

        // For a 16-bit magnitude in a 32-bit container:
        // hr = clz(ax) - (32 - 16) - 1 = clz(ax) - 17
        int hr = static_cast<int>(__builtin_clz(ax)) - 17;
        if (hr < 0) hr = 0;
        if (hr > 14) hr = 14;
        return static_cast<int8_t>(hr);
    }

    /**
     * @brief   Update the saved post-shift value.
     *
     * @param new_shift The new post-shift value.
     */
    void update_post_shift(int8_t new_shift) {
        if (new_shift > _post_shift) {
            const int8_t delta_shift = new_shift - _post_shift;
            _post_shift = new_shift;
            const uint32_t already_set_block_size = _num_biquad_blocks_set * coeffs_per_stage<q15_t>::value;
            arm_shift_q15(_coefficients, static_cast<int8_t>(-delta_shift),
                          _coefficients, already_set_block_size);
        }
    }

    arm_biquad_casd_df1_inst_q15 _cascade_instance{};
    q15_t *_coefficients = nullptr;
    uint32_t _num_biquad_blocks_set = 0;
    int8_t _post_shift = 0;
    int8_t _gain_shift = 0;
};

} // namespace tiny_iir