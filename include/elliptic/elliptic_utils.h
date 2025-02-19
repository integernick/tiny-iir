#pragma once

#include <common/utils.h>

namespace tiny_iir {

static constexpr uint32_t NUM_OF_LANDEN_ITERATIONS = 5;

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
    constexpr double EPS = 1e-6;
    if (k < EPS) {
        k = EPS;
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

/**
 * @brief   Initialize the Landen sequence.
 *
 * @param seq Pointer to the sequence array.
 * @param k   Elliptic modulus.
 */
static void init_landen_sequence(double *seq, double k) {
    seq[0] = landen_next(k);
    for (int i = 1; i < NUM_OF_LANDEN_ITERATIONS; ++i) {
        seq[i] = landen_next(seq[i - 1]);
    }
}

/**
 * @brief   Calculate sn elliptic Jacobi function sn(u, k).
 *
 * @param u Real argument.
 * @param k Elliptic modulus.
 * @return  The elliptic Jacobi function sn(u, k) value.
 */
[[nodiscard]] static Complex sn(Complex u, double k) {
    Complex w = std::sin(u * M_PI_2);

    double landen_seq[NUM_OF_LANDEN_ITERATIONS];
    init_landen_sequence(landen_seq, k);
    for (int i = NUM_OF_LANDEN_ITERATIONS - 1; i >= 0; --i) {
        const double &v_n = landen_seq[i];
        w = (1.0 + v_n) * w / (1.0 + v_n * w * w);
    }

    return w;
}

/**
 * @brief   Calculate sn elliptic Jacobi function cd(u, k).
 *
 * @param u Real argument.
 * @param k Elliptic modulus.
 * @return  The elliptic Jacobi function cd(u, k) value.
 */
[[nodiscard]] Complex cd(Complex u, double k) {
    return sn(u + 1.0, k); // Using cd(z, k) = sn(z + K, k)
}

/**
 * @brief   Solve the degree equation for the elliptic modulus.
 *
 * @param N Filter order.
 * @param k1_prime Complement of (eps_p / eps_s).
 * @return  The ratio of (wp / ws).
 */
[[nodiscard]] double solve_degree_equation(uint32_t N, double k1_prime) {
    const uint32_t L = N / 2;
    double k_prime = std::pow(k1_prime, N);

    for (int i = 0; i < L; ++i) {
        const double u_i = (2.0 * i + 1.0) / N;
        const double sn_val = sn(u_i, k1_prime).real();
        k_prime *= sn_val;
    }

    k_prime = std::pow(k_prime, 4);
    return get_complimentary(k_prime);
}

/**
 * @brief   Calculate inverse of the elliptic Jacobi function sn (complex version).
 *
 * @param w Complex value of the elliptic Jacobi function sn(u, k).
 * @param k Elliptic modulus.
 * @param R The ratio K_prime / K, where K is the complete elliptic integral of the elliptic modulus
 *          and K_prime is the elliptic integral of complement of the elliptic modulus.
 * @return  The argument u of the elliptic Jacobi function w = sn(u, k) value.
 */
[[nodiscard]] Complex asn(Complex w, double k, double R) {
    double v_prev;  // v_{n-1}
    double v_n = k; // v_{n}

    for (int i = 0; i < NUM_OF_LANDEN_ITERATIONS; ++i) {
        v_prev = v_n;
        v_n = landen_next(v_n);
        w = w / (1.0 + std::sqrt(1.0 - w * w * v_prev * v_prev)) * 2.0 / (1.0 + v_n);
    }

    const Complex u = w == 1.0 ? 0.0 : M_2_PI * std::acos(w);

    const Complex acd = {srem(u.real(), 4), srem(u.imag(), 2 * R)};
    return 1.0 - acd; // Using cd(z, k) = sn(z + K, k)
}

}