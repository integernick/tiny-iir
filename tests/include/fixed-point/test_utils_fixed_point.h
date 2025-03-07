#pragma once

#include <type_utils.h>

#include <gtest/gtest.h>

namespace tiny_iir {

static constexpr double TOL_Q31 = 2e-8;

inline void normalize_coeffs(std::vector<double> &coeffs) {
    constexpr uint32_t COEFFICIENTS_PER_BIQUAD_BLOCK = 6;
    const uint32_t num_of_biquad_blocks = coeffs.size() / COEFFICIENTS_PER_BIQUAD_BLOCK;

    for (uint32_t i = 0; i < num_of_biquad_blocks; ++i) {
        std::span<double> biquad_coeffs(coeffs.begin() + i * COEFFICIENTS_PER_BIQUAD_BLOCK,
                                        COEFFICIENTS_PER_BIQUAD_BLOCK);

        const double abs_max_coeff = std::abs(*std::max_element(biquad_coeffs.begin(),
                                                                biquad_coeffs.end(),
                                                                [](double a, double b) {
                                                                    return std::abs(a) < std::abs(b);
                                                                }));
        if (abs_max_coeff > 1.0) {
            for (auto &coeff: biquad_coeffs) {
                coeff /= abs_max_coeff;
            }
        }
    }
}

}