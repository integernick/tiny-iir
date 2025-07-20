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
class BiquadCascade {
public:
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
    void init(uint8_t num_stages,
                     T *coefficients,
                     T *delay_pipeline) {
        biquad_cascade_init(&_biquad_block, num_stages, coefficients, delay_pipeline);
    }

    /**
     * @brief   Process a batch of samples.
     *
     * @param src  The input buffer.
     * @param dst  The output buffer.
     * @param block_size  The number of samples to process.
     */
    void process_cascade(T *src, T *dst, uint32_t block_size) {
        biquad_cascade_process(&_biquad_block, src, dst, block_size);
    }

    /**
     * @brief   Push biquad coefficients to the filter.
     *
     * @param coefficients  The pointer to the biquad coefficients block.
     * @param biquad_coefficients  The biquad coefficients.
     */
    void push_biquad_coefficients(T *coefficients, const BiquadCoefficients &biquad_coefficients) {
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
            //normalized_biquad_coefficients.a1 *= scale;
            //normalized_biquad_coefficients.a2 *= scale;

            /* Accumulate the post shift for the next stage */
            _normalizing_gain *= max_val;
            std::cout << "max_val = " << max_val << "\r\n";
            std::cout << "scale = " << scale << "\r\n";
            std::cout << "normalizing_gain = " << _normalizing_gain << "\r\n";
        }

        coefficients[0] = static_cast<T>(normalized_biquad_coefficients.b0);
        coefficients[1] = static_cast<T>(normalized_biquad_coefficients.b1);
        coefficients[2] = static_cast<T>(normalized_biquad_coefficients.b2);
        coefficients[3] = static_cast<T>(normalized_biquad_coefficients.a1);
        coefficients[4] = static_cast<T>(normalized_biquad_coefficients.a2);
        std::cout << "push biquad coefficients\r\n";
        std::cout << "b0 = " << coefficients[0] << "\r\n";
        std::cout << "b1 = " << coefficients[1] << "\r\n";
        std::cout << "b2 = " << coefficients[2] << "\r\n";
        std::cout << "a1 = " << coefficients[3] << "\r\n";
        std::cout << "a2 = " << coefficients[4] << "\r\n";
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = BiquadBlockDF2TCascade<T>::BLOCK_DELAY_LINE_SIZE;

private:
    BiquadBlockDF2TCascade<T> _biquad_block;
    double _normalizing_gain = 1.0;
};

}
