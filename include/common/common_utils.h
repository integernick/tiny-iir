#pragma once

#include <complex>

namespace tiny_iir {

using Complex = std::complex<double>;

static constexpr double INFINITY_VALUE
        = std::numeric_limits<double>::infinity();

enum class FilterPassType : int {
    LOW_PASS,
    HIGH_PASS,
    BAND_PASS,
    BAND_STOP,
};

struct PoleZeroPair {
    std::complex<double> pole;
    std::complex<double> zero;
};

struct BiquadCoefficients {
    double b0 = 0;
    double b1 = 0;
    double b2 = 0;
    double a1 = 0;
    double a2 = 0;
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
T constrain(T value, T min, T max) {
    if (value < min) {
        return min;
    } else if (value > max) {
        return max;
    } else {
        return value;
    }
}

}
