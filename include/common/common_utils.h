#pragma once

#include <complex>

namespace tiny_iir {

template<typename DESIGN_T>
using Complex = std::complex<DESIGN_T>;

enum class FilterPassType : int {
    LowPass,
    HighPass,
    BandPass,
    BandStop,

    LOW_PASS = LowPass,
    HIGH_PASS = HighPass,
    BAND_PASS = BandPass,
    BAND_STOP = BandStop
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
 * @tparam T The type of the value to constrain.
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

/**
 * @brief   Clamp frequency to safe range..
 *
 * @details Frequencies too close to 0 or Nyquist can cause numerical issues
 *          in the bilinear transform.
 *
 * @tparam DESIGN_T         The design type (float or double).
 * @param normalized_freq   Normalized frequency (0 to 1, where 1 = Nyquist).
 * @return  Clamped frequency in safe range.
 */
template<typename DESIGN_T>
DESIGN_T clamp_frequency(DESIGN_T frequency) {
    constexpr DESIGN_T MIN_FREQ = DESIGN_T{0.0001};  // 0.01% of Nyquist
    constexpr DESIGN_T MAX_FREQ = DESIGN_T{0.9999};  // 99.99% of Nyquist

    return constrain(std::abs(frequency), MIN_FREQ, MAX_FREQ);
}

}
