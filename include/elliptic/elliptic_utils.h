#pragma once

#include <common/common_utils.h>
#include <numbers>

namespace tiny_iir {

static constexpr uint32_t NUM_OF_LANDEN_ITERATIONS = 5;

/**
 * @brief   Get complimentary of elliptic modulus.
 *
 * @tparam  DESIGN_T  The filter design type.
 * @param k Elliptic modulus.
 * @return  The complimentary of the elliptic modulus.
 */
template<typename DESIGN_T>
[[nodiscard]] inline DESIGN_T get_complimentary(DESIGN_T k) {
    return std::sqrt(DESIGN_T{1} - k * k);
}

/**
 * @brief   Calculate arithmetic-geometric mean.
 *
 * @tparam  DESIGN_T  The filter design type.
 * @param a First value.
 * @param b Second value.
 * @return  Arithmetic-geometric mean of a and b.
 */
template<typename DESIGN_T>
[[nodiscard]] inline DESIGN_T arithmetic_geometric_mean(DESIGN_T a, DESIGN_T b) {
    if (a < b) {
        std::swap(a, b);
    }

    DESIGN_T c = a - b;
    DESIGN_T c_prev;

    do {
        c_prev = c;
        c = (a - b) / DESIGN_T{2};
        DESIGN_T a_new = (a + b) / DESIGN_T{2};
        b = std::sqrt(a * b);
        a = a_new;
    } while (c < c_prev);

    return a;
}

/**
 * @brief   Calculate Landen transformation at step n+1.
 *
 * @tparam  DESIGN_T  The filter design type.
 * @param k Elliptic modulus at step n.
 * @return  Elliptic modulus at step n+1.
 */
template<typename DESIGN_T>
[[nodiscard]] inline DESIGN_T landen_next(DESIGN_T k) {
    const DESIGN_T k_comp = get_complimentary(k);
    return (DESIGN_T{1} - k_comp) / (DESIGN_T{1} + k_comp);
}

/**
 * @brief   Calculate elliptic integral of the first kind.
 *
 * @tparam  DESIGN_T  The filter design type.
 * @param k1    Elliptic modulus.
 * @return  Value of the elliptic integral of the first kind.
 */
template<typename DESIGN_T>
[[nodiscard]] inline DESIGN_T calculate_elliptic_integral(DESIGN_T k) {
    constexpr auto EPS = static_cast<DESIGN_T>(1e-6);

    if (k < EPS) {
        k = EPS;
    }

    const DESIGN_T agm = arithmetic_geometric_mean(DESIGN_T{1}, get_complimentary(k));
    return std::numbers::pi_v<DESIGN_T> / (DESIGN_T{2} * agm);
}

/**
 * @brief   Symmetric remainder. Returns a value z such that z is in [-y/2, y/2].
 *
 * @tparam  DESIGN_T  The filter design type.
 * @param x A real-valued floating point number.
 * @param y A positive floating point number.
 * @return The symmetric remainder of x w.r.t. y.
 */
template<typename DESIGN_T>
[[nodiscard]] inline DESIGN_T srem(DESIGN_T x, DESIGN_T y) {
    // First, bring x into the interval (-y, y] using fmod:
    DESIGN_T z = std::fmod(x, y);

    // Now, if z is outside [-y/2, y/2], shift it by ±y:
    if (z > y / DESIGN_T{2}) {
        z -= y;
    } else if (z <= -y / DESIGN_T{2}) {
        z += y;
    }

    return z;
}

/**
 * @brief   Initialize the Landen sequence.
 *
 * @tparam  DESIGN_T  The filter design type.
 * @param seq Pointer to the sequence array.
 * @param k   Elliptic modulus.
 */
template<typename DESIGN_T>
static inline void init_landen_sequence(DESIGN_T *seq, DESIGN_T k) {
    seq[0] = landen_next(k);

    for (uint32_t i = 1; i < NUM_OF_LANDEN_ITERATIONS; ++i) {
        seq[i] = landen_next(seq[i - 1]);
    }
}

/**
 * @brief   Calculate sn elliptic Jacobi function sn(u, k).
 *
 * @tparam  DESIGN_T  The filter design type.
 * @param u Real argument.
 * @param k Elliptic modulus.
 * @return  The elliptic Jacobi function sn(u, k) value.
 */
template<typename DESIGN_T>
[[nodiscard]] inline static Complex<DESIGN_T> sn(Complex<DESIGN_T> u, DESIGN_T k) {
    Complex<DESIGN_T> w = std::sin(u * std::numbers::pi_v<DESIGN_T> / 2.0);

    DESIGN_T landen_seq[NUM_OF_LANDEN_ITERATIONS];
    init_landen_sequence(landen_seq, k);

    for (int i = NUM_OF_LANDEN_ITERATIONS - 1; i >= 0; --i) {
        const DESIGN_T &v_n = landen_seq[i];
        w = (DESIGN_T{1} + v_n) * w / (DESIGN_T{1} + v_n * w * w);
    }

    return w;
}

/**
 * @brief   Calculate sn elliptic Jacobi function cd(u, k).
 *
 * @tparam  DESIGN_T  The filter design type.
 * @param u Real argument.
 * @param k Elliptic modulus.
 * @return  The elliptic Jacobi function cd(u, k) value.
 */
template<typename DESIGN_T>
[[nodiscard]] inline Complex<DESIGN_T> cd(Complex<DESIGN_T> u, DESIGN_T k) {
    return sn<DESIGN_T>(u + DESIGN_T{1}, k); // Using cd(z, k) = sn(z + K, k)
}

/**
 * @brief   Solve the degree equation for the elliptic modulus.
 *
 * @tparam  DESIGN_T  The filter design type.
 * @param N Filter order.
 * @param k1_prime Complement of (eps_p / eps_s).
 * @return  The ratio of (wp / ws).
 */
template<typename DESIGN_T>
[[nodiscard]] inline DESIGN_T solve_degree_equation(uint32_t N, DESIGN_T k1_prime) {
    const uint32_t L = N / 2;
    DESIGN_T k_prime = std::pow(k1_prime, N);

    for (uint32_t i = 0; i < L; ++i) {
        const DESIGN_T u_i = (DESIGN_T{2} * i + DESIGN_T{1}) / N;
        const DESIGN_T sn_val = sn<DESIGN_T>(u_i, k1_prime).real();
        k_prime *= sn_val;
    }

    k_prime = std::pow(k_prime, 4);
    return get_complimentary(k_prime);
}

/**
 * @brief   Calculate inverse of the elliptic Jacobi function sn (complex version).
 *
 * @tparam  DESIGN_T  The filter design type.
 * @param w Complex value of the elliptic Jacobi function sn(u, k).
 * @param k Elliptic modulus.
 * @param R The ratio K_prime / K, where K is the complete elliptic integral of the elliptic modulus
 *          and K_prime is the elliptic integral of complement of the elliptic modulus.
 * @return  The argument u of the elliptic Jacobi function w = sn(u, k) value.
 */
template<typename DESIGN_T>
[[nodiscard]] inline Complex<DESIGN_T> asn(Complex<DESIGN_T> w, DESIGN_T k, DESIGN_T R) {
    DESIGN_T v_prev;  // v_{n-1}
    DESIGN_T v_n = k; // v_{n}

    for (uint32_t i = 0; i < NUM_OF_LANDEN_ITERATIONS; ++i) {
        v_prev = v_n;
        v_n = landen_next(v_n);
        w = w / (DESIGN_T{1} + std::sqrt(DESIGN_T{1} - w * w * v_prev * v_prev)) * DESIGN_T{2} / (DESIGN_T{1} + v_n);
    }

    const Complex u = (w == DESIGN_T{1})
                      ? DESIGN_T{0}
                      : DESIGN_T{2} * std::numbers::inv_pi * std::acos(w);

    const Complex acd = {srem<DESIGN_T>(u.real(), DESIGN_T{4}), srem(u.imag(), DESIGN_T{2} * R)};
    return DESIGN_T{1} - acd; // Using cd(z, k) = sn(z + K, k)
}

}