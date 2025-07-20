#pragma once

#include <type_utils.h>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

namespace tiny_iir {

constexpr double TOL_Q31 = 2e-8;

inline void normalize_coeffs(std::vector<double> &coeffs) {
    constexpr uint32_t COEFFS_PER_BLOCK = 6; // b0,b1,b2,a0,a1,a2
    const uint32_t blocks = coeffs.size() / COEFFS_PER_BLOCK;

    // First pass: global max across all stages (exclude a0 from max if you prefer; including it is harmless)
    double max_abs = 0.0;
    for (uint32_t i = 0; i < blocks; ++i) {
        const double *blk = coeffs.data() + i * COEFFS_PER_BLOCK;
        const double b0 = blk[0];
        const double b1 = blk[1];
        const double b2 = blk[2];
        const double a1 = blk[4];
        const double a2 = blk[5];
        max_abs = std::max({max_abs, std::abs(b0), std::abs(b1), std::abs(b2), std::abs(a1), std::abs(a2)});
    }

    // Choose k so that max_abs / 2^k < 1.0  => k = ceil(log2(max_abs)) when max_abs>1
    int k = 0;
    if (max_abs > 1.0) {
        k = static_cast<int>(std::ceil(std::log2(max_abs)));
        if (k < 0) k = 0;
    }
    const double scale = std::ldexp(1.0, -k); // 2^{-k}

    // Second pass: scale all six entries per block by 2^{-k}
    for (uint32_t i = 0; i < blocks; ++i) {
        double *blk = coeffs.data() + i * COEFFS_PER_BLOCK;
        blk[0] *= scale; // b0
        blk[1] *= scale; // b1
        blk[2] *= scale; // b2
        blk[3] *= scale; // a0 (harmless to scale for expected values)
        blk[4] *= scale; // a1
        blk[5] *= scale; // a2
    }
}

}