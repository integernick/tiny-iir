#pragma once

#include <complex>

namespace tiny_iir {

template<typename DESIGN_T>
using Complex = std::complex<DESIGN_T>;

enum class FilterPassType : int {
    LOW_PASS,
    HIGH_PASS,
    BAND_PASS,
    BAND_STOP,
};

template<typename DESIGN_T>
struct PoleZeroPair {
    std::complex<DESIGN_T> pole;
    std::complex<DESIGN_T> zero;
};

template<typename DESIGN_T>
struct BiquadCoefficients {
    DESIGN_T b0 = 0;
    DESIGN_T b1 = 0;
    DESIGN_T b2 = 0;
    DESIGN_T a1 = 0;
    DESIGN_T a2 = 0;
};

/**
 * @brief   Constrain a value to a specified range.
 *
 * @param value  The value to constrain.
 * @param min  The minimum value.
 * @param max  The maximum value.
 * @return  The constrained value.
 */
template <typename T>
constexpr T constrain(T value, T min, T max) {
    return value < min ? min
            : value > max ? max
            : value;
}

}
