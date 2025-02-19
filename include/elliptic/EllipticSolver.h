#pragma once

#include <common/utils.h>
#include "elliptic_utils.h"

namespace tiny_iir {

class EllipticSolver {
public:
    /**
     * @brief   Initialize the elliptic solver: calculate complete elliptic integrals.
     *
     * @param k Elliptic modulus.
     */
    void init(double k) {
        _k = k;
        _k_prime = std::sqrt(1.0 - _k * _k);
        _K = calculate_elliptic_integral(_k);
        _K_prime = calculate_elliptic_integral(_k_prime);
    }

    /**
     * @brief   Calculate sn elliptic Jacobi function.
     *
     * @param u Argument.
     * @return  Elliptic modulus at step n.
     */
    [[nodiscard]] double sn(double u) const {
        return sn(u, _k);
    }

    /**
     * @brief   Calculate cd elliptic Jacobi function (complex version).
     *
     * @param u Complex argument.
     * @return  The elliptic Jacobi function cd value.
     */
    [[nodiscard]] Complex cd(Complex u) const {
        Complex w = std::cos(u * M_PI_2);

        double v_n = _k;
        for (int i = 0; i < NUM_OF_LANDEN_ITERATIONS; ++i) {
            v_n = landen_next(v_n);
            w = (1.0 + v_n) * w / (1.0 + v_n * w * w);
        }

        return w;
    }

    /**
     * @brief   Calculate cd elliptic Jacobi function (real version).
     *
     * @param u Real argument.
     * @return  The elliptic Jacobi function cd value.
     */
    [[nodiscard]] double cd(double u) const {
        double w = std::cos(u * M_PI_2);

        double v_n = _k;
        for (int i = 0; i < NUM_OF_LANDEN_ITERATIONS; ++i) {
            v_n = landen_next(v_n);
            w = (1.0 + v_n) * w / (1.0 + v_n * w * w);
        }

        return w;
    }

    /**
     * @brief   Calculate inverse of the elliptic Jacobi function cd (complex version).
     *
     * @param w Complex value of the elliptic Jacobi function cd(u, k).
     * @return  The argument u of the elliptic Jacobi function w = cd(u, k) value.
     */
    [[nodiscard]] Complex acd(Complex w) const {
        double v_prev;  // v_{n-1}
        double v_n = _k; // v_{n}

        for (int i = 0; i < NUM_OF_LANDEN_ITERATIONS; ++i) {
            v_prev = v_n;
            v_n = landen_next(v_n);
            w = w / (1.0 + std::sqrt(1.0 - w * w * v_prev * v_prev)) * 2.0 / (1.0 + v_n);
        }

        const Complex u = w == 1.0 ? 0.0 : M_2_PI * std::acos(w);
        const double R = _K_prime / _K;

        return {srem(u.real(), 4), srem(u.imag(), 2 * R)};
    }

    /**
     * @brief   Calculate inverse of the elliptic Jacobi function cd (real version).
     *
     * @param w Real value of the elliptic Jacobi function cd(u, k).
     * @return  The argument u of the elliptic Jacobi function w = cd(u, k) value.
     */
    [[nodiscard]] double acd(double w) const {
        double v_prev;  // v_{n-1}
        double v_n = _k; // v_{n}

        for (int i = 0; i < NUM_OF_LANDEN_ITERATIONS; ++i) {
            v_prev = v_n;
            v_n = landen_next(v_n);
            w = w / (1.0 + std::sqrt(1.0 - w * w * v_prev * v_prev)) * 2.0 / (1.0 + v_n);
        }

        const double u = w == 1.0 ? 0.0 : M_2_PI * std::acos(w);

        return u;
    }

    /**
     * @brief   Calculate inverse of the elliptic Jacobi function sn (complex version).
     *
     * @param w Complex value of the elliptic Jacobi function sn(u, k).
     * @return  The argument u of the elliptic Jacobi function w = sn(u, k) value.
     */
    [[nodiscard]] Complex asn(Complex w) const {
        return 1.0 - acd(w);
    }

    /**
     * @brief   Solve the degree equation for the elliptic modulus.
     *
     * @param N Number of sections.
     * @return  The elliptic modulus.
     */
    [[nodiscard]] double asn(double w) const {
        return 1.0 - acd(w);
    }

    /**
     * @brief   Solve the degree equation for the elliptic modulus.
     *
     * @param N Filter order.
     * @return  The elliptic modulus.
     */
    [[nodiscard]] double solve_degree_equation(uint32_t N) const {
        const uint32_t L = N / 2;
        double k_comp = std::pow(_k_prime, N);

        for (int i = 0; i < L; ++i) {
            const double u_i = (2.0 * i + 1.0) / N;
            const double sn_val = sn(u_i, _k_prime);
            k_comp *= sn_val;
        }

        k_comp = std::pow(k_comp, 4);
        return get_complimentary(k_comp);
    }

private:
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
     * @brief   Calculate sn elliptic Jacobi function (real version).
     *
     * @param u Real argument.
     * @param k Elliptic modulus.
     * @return  The elliptic Jacobi function sn(u, k) value.
     */
    [[nodiscard]] static double sn(double u, double k) {
        double w = std::sin(u * M_PI_2);

        double landen_seq[NUM_OF_LANDEN_ITERATIONS];
        init_landen_sequence(landen_seq, k);
        for (int i = NUM_OF_LANDEN_ITERATIONS - 1; i >= 0; --i) {
            const double &v_n = landen_seq[i];
            w = (1.0 + v_n) * w / (1.0 + v_n * w * w);
        }

        return w;
    }

    static constexpr uint32_t NUM_OF_LANDEN_ITERATIONS = 5;

    double _k = 0;
    double _k_prime = 0;
    double _K = 0;
    double _K_prime = 0;
};

}