#pragma once

#include "BiquadBlock.h"

#include <common/common_utils.h>

namespace tiny_iir {

using q31_t = int32_t;
using q15_t = int16_t;

/**
 * @brief   Biquad cascade for double native type.
 */
template<typename T>
struct coeffs_per_stage {
    static constexpr uint32_t value = 5;
};

template<typename T, typename DESIGN_T>
struct BiquadCascade;

/**
 * @brief   Biquad cascade for float native type.
 */
template<typename T, typename DESIGN_T>
class BiquadCascade {
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
    void set_coefficients(T *coefficients) {
        _coefficients = coefficients;
    }

    /**
     * @brief   Initialize the biquad cascade filter.
     *
     * @param num_stages  The number of stages.
     * @param delay_pipeline  The delay pipeline array.
     */
    void init(uint8_t num_stages, T *delay_pipeline) {
        biquad_cascade_init(&_biquad_cascade_instance, num_stages, _coefficients, delay_pipeline);
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
    void process_cascade(T *src, T *dst, uint32_t block_size) {
        biquad_cascade_process(&_biquad_cascade_instance, src, dst, block_size);
    }

    /**
     * @brief   Push biquad coefficients to the filter.
     *
     * @param coefficients  The pointer to the biquad coefficients block.
     * @param biquad_coefficients  The biquad coefficients.
     */
    void push_biquad_coefficients(T *coefficients, const BiquadCoefficients<DESIGN_T> &biquad_coefficients) {
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
    [[nodiscard]] BiquadCoefficients<DESIGN_T> get_biquad_coefficients(uint32_t biquad_idx) const {
        BiquadCoefficients<DESIGN_T> coeffs;
        to_native<DESIGN_T>(&_coefficients[biquad_idx * coeffs_per_stage<float>::value], &coeffs.b0,
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
    void update_gain(DESIGN_T biquad_gain_inv, T &gain) {
        gain *= static_cast<T>(biquad_gain_inv);
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = BiquadBlockDF2TCascade<T>::BLOCK_DELAY_LINE_SIZE;

    static constexpr double UNITY = T{1};

private:
    BiquadBlockDF2TCascade<T> _biquad_cascade_instance;
    T *_coefficients = nullptr;
    uint32_t _num_biquad_blocks_set = 0;
};

}
