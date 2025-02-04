#pragma once

#include <complex>

namespace tiny_iir {

using Complex = std::complex<double>;

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
