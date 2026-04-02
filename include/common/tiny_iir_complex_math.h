#pragma once

#include "tiny_iir_math.h"

namespace tiny_iir {

// ── Complex transcendentals ─────────────────────────────────────────────────
//
// When TINY_IIR_USE_STL is defined, delegate to the standard <complex> overloads.
// Otherwise, provide lightweight implementations using real-valued <math.h>.

#if defined(TINY_IIR_USE_STL)

template<typename T>
Complex<T> complex_sqrt(const Complex<T> &z) { return std::sqrt(z); }

template<typename T>
Complex<T> complex_sin(const Complex<T> &z) { return std::sin(z); }

template<typename T>
Complex<T> complex_log(const Complex<T> &z) { return std::log(z); }

template<typename T>
Complex<T> complex_acos(const Complex<T> &z) { return std::acos(z); }

#else

/**
 * @brief   Complex square root.
 *
 * @details Uses the standard decomposition into real-valued operations.
 *          Branch cut along the negative real axis.
 */
template<typename T>
Complex<T> complex_sqrt(const Complex<T> &z) {
    const T a = z.re;
    const T b = z.im;

    const T r = abs(z);

    if (r == T{0}) {
        return {T{0}, T{0}};
    }

    if (a >= T{0}) {
        const T t = static_cast<T>(sqrt(static_cast<double>(r + a) * 0.5));
        return {t, b / (T{2} * t)};
    } else {
        const T t = static_cast<T>(sqrt(static_cast<double>(r - a) * 0.5));
        return {static_cast<T>(fabs(static_cast<double>(b))) / (T{2} * t),
                static_cast<T>(copysign(static_cast<double>(t), static_cast<double>(b)))};
    }
}

/**
 * @brief   Complex sine.
 *
 * @details sin(a + bi) = sin(a)cosh(b) + i*cos(a)sinh(b)
 */
template<typename T>
Complex<T> complex_sin(const Complex<T> &z) {
    const auto a = static_cast<double>(z.re);
    const auto b = static_cast<double>(z.im);
    return {static_cast<T>(sin(a) * cosh(b)),
            static_cast<T>(cos(a) * sinh(b))};
}

/**
 * @brief   Complex natural logarithm.
 *
 * @details log(z) = 0.5 * log(|z|^2) + i * atan2(im, re)
 */
template<typename T>
Complex<T> complex_log(const Complex<T> &z) {
    const auto a = static_cast<double>(z.re);
    const auto b = static_cast<double>(z.im);
    return {static_cast<T>(0.5 * log(a * a + b * b)),
            static_cast<T>(atan2(b, a))};
}

/**
 * @brief   Complex arc cosine.
 *
 * @details acos(z) = -i * log(z + i * sqrt(1 - z^2))
 */
template<typename T>
Complex<T> complex_acos(const Complex<T> &z) {
    const Complex<T> i_unit = {T{0}, T{1}};
    const Complex<T> one_minus_z2 = Complex<T>{T{1}} - z * z;
    const Complex<T> s = complex_sqrt(one_minus_z2);
    const Complex<T> arg = z + i_unit * s;
    const Complex<T> lg = complex_log(arg);
    // -i * lg = {lg.im, -lg.re}
    return {lg.im, -lg.re};
}

#endif // !TINY_IIR_USE_STL

} // namespace tiny_iir
