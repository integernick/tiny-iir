#pragma once

/**
 * @file    tiny_iir_math.h
 * @brief   Math primitives for tiny-iir.
 *
 * @details By default, provides lightweight freestanding implementations of
 *          Complex<T>, math constants, and utilities — no hosted C++ headers needed.
 *
 *          Define TINY_IIR_USE_STL (or pass -DTINY_IIR_USE_STL) to use the
 *          standard library (<complex>, <numbers>, <algorithm>, <limits>) instead.
 */

#include <math.h>
#include <type_traits>

#if defined(TINY_IIR_USE_STL)
#include <algorithm>
#include <complex>
#include <limits>
#include <numbers>
#include <utility>
#endif

namespace tiny_iir {

// ── Math constants ──────────────────────────────────────────────────────────

namespace numbers {

#if defined(TINY_IIR_USE_STL)

template<typename T>
inline constexpr T pi_v = std::numbers::pi_v<T>;

template<typename T>
inline constexpr T ln10_v = std::numbers::ln10_v<T>;

template<typename T>
inline constexpr T inv_pi_v = std::numbers::inv_pi_v<T>;

#else

template<typename T>
inline constexpr T pi_v = T(3.14159265358979323846L);

template<typename T>
inline constexpr T ln10_v = T(2.30258509299404568402L);

template<typename T>
inline constexpr T inv_pi_v = T(0.31830988618379067154L);

#endif

} // namespace numbers

// ── Infinity ────────────────────────────────────────────────────────────────

#if defined(TINY_IIR_USE_STL)

template<typename T>
inline constexpr T infinity_v = std::numeric_limits<T>::infinity();

#else

template<typename T>
inline constexpr T infinity_v = T(INFINITY);

#endif

// ── Complex<T> ──────────────────────────────────────────────────────────────

#if defined(TINY_IIR_USE_STL)

template<typename T>
using Complex = std::complex<T>;

using std::abs;
using std::norm;
using std::conj;
using std::polar;

#else

template<typename T>
struct Complex {
    T re{};
    T im{};

    constexpr Complex() = default;
    constexpr Complex(T real, T imag) : re(real), im(imag) {}

    // Implicit scalar-to-complex conversion (critical for expressions like DT{1} + s)
    constexpr Complex(T real) : re(real), im(T{0}) {} // NOLINT(google-explicit-constructor)

    [[nodiscard]] constexpr T real() const { return re; }
    [[nodiscard]] constexpr T imag() const { return im; }

    constexpr void real(T r) { re = r; }
    constexpr void imag(T i) { im = i; }

    // Compound assignment — complex
    constexpr Complex &operator+=(const Complex &rhs) {
        re += rhs.re;
        im += rhs.im;
        return *this;
    }

    constexpr Complex &operator-=(const Complex &rhs) {
        re -= rhs.re;
        im -= rhs.im;
        return *this;
    }

    constexpr Complex &operator*=(const Complex &rhs) {
        const T r = re * rhs.re - im * rhs.im;
        const T i = re * rhs.im + im * rhs.re;
        re = r;
        im = i;
        return *this;
    }

    constexpr Complex &operator/=(const Complex &rhs) {
        const T denom = rhs.re * rhs.re + rhs.im * rhs.im;
        const T r = (re * rhs.re + im * rhs.im) / denom;
        const T i = (im * rhs.re - re * rhs.im) / denom;
        re = r;
        im = i;
        return *this;
    }

    // Compound assignment — scalar
    constexpr Complex &operator*=(T s) {
        re *= s;
        im *= s;
        return *this;
    }

    constexpr Complex &operator/=(T s) {
        re /= s;
        im /= s;
        return *this;
    }

    // Unary minus
    constexpr Complex operator-() const {
        return {-re, -im};
    }

    // Comparison
    constexpr bool operator==(const Complex &rhs) const {
        return re == rhs.re && im == rhs.im;
    }

    constexpr bool operator!=(const Complex &rhs) const {
        return !(*this == rhs);
    }
};

// ── Complex binary operators ────────────────────────────────────────────────

// Complex OP Complex
template<typename T>
constexpr Complex<T> operator+(Complex<T> a, const Complex<T> &b) { a += b; return a; }

template<typename T>
constexpr Complex<T> operator-(Complex<T> a, const Complex<T> &b) { a -= b; return a; }

template<typename T>
constexpr Complex<T> operator*(Complex<T> a, const Complex<T> &b) { a *= b; return a; }

template<typename T>
constexpr Complex<T> operator/(Complex<T> a, const Complex<T> &b) { a /= b; return a; }

// Scalar OP Complex
template<typename T>
constexpr Complex<T> operator+(T s, const Complex<T> &c) { return Complex<T>(s) + c; }

template<typename T>
constexpr Complex<T> operator-(T s, const Complex<T> &c) { return Complex<T>(s) - c; }

template<typename T>
constexpr Complex<T> operator*(T s, const Complex<T> &c) { return {s * c.re, s * c.im}; }

template<typename T>
constexpr Complex<T> operator/(T s, const Complex<T> &c) { return Complex<T>(s) / c; }

// Complex OP Scalar
template<typename T>
constexpr Complex<T> operator+(const Complex<T> &c, T s) { return {c.re + s, c.im}; }

template<typename T>
constexpr Complex<T> operator-(const Complex<T> &c, T s) { return {c.re - s, c.im}; }

template<typename T>
constexpr Complex<T> operator*(const Complex<T> &c, T s) { return {c.re * s, c.im * s}; }

template<typename T>
constexpr Complex<T> operator/(const Complex<T> &c, T s) { return {c.re / s, c.im / s}; }

// ── Complex free functions ──────────────────────────────────────────────────

template<typename T>
T abs(const Complex<T> &z) {
    return static_cast<T>(sqrt(static_cast<double>(z.re) * z.re + static_cast<double>(z.im) * z.im));
}

template<typename T>
constexpr T norm(const Complex<T> &z) {
    return z.re * z.re + z.im * z.im;
}

template<typename T>
constexpr Complex<T> conj(const Complex<T> &z) {
    return {z.re, -z.im};
}

template<typename T>
Complex<T> polar(T mag, T angle) {
    return {mag * static_cast<T>(cos(static_cast<double>(angle))),
            mag * static_cast<T>(sin(static_cast<double>(angle)))};
}

#endif // !TINY_IIR_USE_STL

// ── Utility functions ───────────────────────────────────────────────────────

#if defined(TINY_IIR_USE_STL)

template<typename T>
constexpr T tiny_iir_min(T a, T b) { return std::min(a, b); }

template<typename T>
constexpr T tiny_iir_max(T a, T b) { return std::max(a, b); }

template<typename T, typename... Args>
constexpr T tiny_iir_max(T a, T b, Args... rest) { return tiny_iir_max(std::max(a, b), rest...); }

template<typename T>
constexpr void tiny_iir_swap(T &a, T &b) { std::swap(a, b); }

template<typename A, typename B>
using Pair = std::pair<A, B>;

#else

template<typename T>
constexpr T tiny_iir_min(T a, T b) {
    return (a < b) ? a : b;
}

template<typename T>
constexpr T tiny_iir_max(T a, T b) {
    return (a > b) ? a : b;
}

template<typename T, typename... Args>
constexpr T tiny_iir_max(T a, T b, Args... rest) { return tiny_iir_max(tiny_iir_max(a, b), rest...); }

template<typename T>
constexpr void tiny_iir_swap(T &a, T &b) {
    T tmp = a;
    a = b;
    b = tmp;
}

template<typename A, typename B>
struct Pair {
    A first;
    B second;
};

#endif // !TINY_IIR_USE_STL

} // namespace tiny_iir
