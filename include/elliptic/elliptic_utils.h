#pragma once

#include <complex>

namespace tiny_iir {

/**
 * @brief   Get complimentary of elliptic modulus.
 *
 * @param k Elliptic modulus.
 * @return  The complimentary of the elliptic modulus.
 */
[[nodiscard]] inline double get_complimentary(double k) {
    return std::sqrt(1.0 - k * k);
}

/**
 * @brief   Calculate arithmetic-geometric mean.
 *
 * @param a First value.
 * @param b Second value.
 * @return  Arithmetic-geometric mean of a and b.
 */
[[nodiscard]] inline double arithmetic_geometric_mean(double a, double b) {
    if (a < b) {
        std::swap(a, b);
    }
    double c = a - b;
    double c_prev;
    do {
        c_prev = c;
        c = (a - b) / 2;
        double a_new = (a + b) / 2;
        b = std::sqrt(a * b);
        a = a_new;
    } while (c < c_prev);

    return a;
}

/**
 * @brief   Calculate Landen transformation at step n+1.
 *
 * @param k Elliptic modulus at step n.
 * @return  Elliptic modulus at step n+1.
 */
[[nodiscard]] inline double landen_next(double k) {
    const double k_comp = get_complimentary(k);
    return (1.0 - k_comp) / (1.0 + k_comp);
}

/**
 * @brief   Calculate elliptic integral of the first kind.
 *
 * @param k1    Elliptic modulus.
 * @return  Value of the elliptic integral of the first kind.
 */
[[nodiscard]] double calculate_elliptic_integral(double k) {
    if (k < 1e-6) {
        k = 1e-6;
    }

    const double agm = arithmetic_geometric_mean(1.0, get_complimentary(k));
    return M_PI / (2 * agm);
}

/**
 * @brief   Symmetric remainder. Returns a value z such that z is in [-y/2, y/2].
 *
 * @param x A real-valued floating point number.
 * @param y A positive floating point number.
 * @return The symmetric remainder of x w.r.t. y.
 */
[[nodiscard]] inline double srem(double x, double y) {
    // First, bring x into the interval (-y, y] using fmod:
    double z = std::fmod(x, y);

    // Now, if z is outside [-y/2, y/2], shift it by ±y:
    if (z > y / 2.0) {
        z -= y;
    } else if (z <= -y / 2.0) {
        z += y;
    }

    return z;
}

}