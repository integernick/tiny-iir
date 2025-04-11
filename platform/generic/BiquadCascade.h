#pragma once

#include <common/common_utils.h>

#include "BiquadBlock.h"

#include <array>

#if CASCADE_FILTER_DEBUG > 0

#include <iostream>
#include <iomanip>

#endif

namespace tiny_iir {

template<typename>
struct BiquadCascade;

/**
 * @brief   Biquad cascade for float native type.
 */
template<typename T = double>
struct BiquadCascade {
    using type = BiquadBlockDF2TCascade<T>;

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
                     T *coefficients,
                     T *delay_pipeline) {
        biquad_cascade_init(cascade, num_stages, coefficients, delay_pipeline);
    }

    /**
     * @brief   Process a batch of samples.
     *
     * @param cascade  The filter instance.
     * @param src  The input buffer.
     * @param dst  The output buffer.
     * @param block_size  The number of samples to process.
     */
    static void process_cascade(type *cascade, T *src, T *dst, uint32_t block_size) {
        biquad_cascade_process(cascade, src, dst, block_size);
    }

    /**
     * @brief   Multiply two numbers (float version).
     *
     * @param x  The first number.
     * @param y  The second number.
     * @return  The product of the two numbers.
     */
    [[nodiscard]] static T multiply(T x, T y) {
        return x * y;
    }

    /**
     * @brief   Push biquad coefficients to the filter.
     *
     * @param coefficients  The pointer to the biquad coefficients block.
     * @param biquad_coefficients  The biquad coefficients.
     */
    static void push_biquad_coefficients(T *coefficients, const BiquadCoefficients &biquad_coefficients) {
        coefficients[0] = static_cast<T>(biquad_coefficients.b0);
        coefficients[1] = static_cast<T>(biquad_coefficients.b1);
        coefficients[2] = static_cast<T>(biquad_coefficients.b2);
        coefficients[3] = static_cast<T>(-biquad_coefficients.a1);
        coefficients[4] = static_cast<T>(-biquad_coefficients.a2);
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = type::BLOCK_DELAY_LINE_SIZE;
};

}
