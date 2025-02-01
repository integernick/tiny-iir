#pragma once

#include <cstring>

namespace tiny_iir {

/**
 * @addtogroup tiny_iir_module
 * @{
 */

/**
 * @brief   Biquad block (direct form I).
 *
 * @details Implements H(z) = Y(z) / X(z) = H1(z) * H2(z), where
 *          H1(z) = 1 / (1 - a1 * z^-1 - a2 * z^-2) = V(z) / X(z),
 *          H2(z) = b0 + b1 * z^-1 + b2 * z^-2 = Y(z) / V(z).
 * @tparam  T       Data type.
 * @tparam  FORM    Block form.
 */
template <typename T = double>
class BiquadBlockDF1 {
public:
    void set_coefficients(const T coefficients[5]) {
        memcpy(_B, coefficients, sizeof(T) * 3);
        memcpy(_A, coefficients + 3, sizeof(T) * 2);
    }

    T process(T input) {
        const T v = input + _A[0] * _delay_line[0] + _A[1] * _delay_line[1];
        const T output = _B[0] * v + _B[1] * _delay_line[0] + _B[2] * _delay_line[1];

        _delay_line[1] = _delay_line[0];
        _delay_line[0] = v;

        return output;
    }

    void reset() {
        _delay_line[0] = _delay_line[1] = 0;
    }

private:
    T _B[3]{}; // b0, b1, b2
    T _A[2]{}; // a1, a2

    T _delay_line[2]{}; // v[n-1], v[n-2]
};

/**
 * @brief   Biquad block (direct form II or canonical form).
 *
 * @details Implements H(z) = (b0 + b1 * z^-1 + b2 * z^-2) / (1 - a1 * z^-1 - a2 * z^-2)
 * @tparam  T       Data type.
 * @tparam  FORM    Block form.
 */
template <typename T = double>
class BiquadBlockDF2 {
public:
    void set_coefficients(const T coefficients[5]) {
        memcpy(_B, coefficients, sizeof(T) * 3);
        memcpy(_A, coefficients + 3, sizeof(T) * 2);
    }

    T process(T input) {
        const T output = _B[0] * input + _B[1] * _delay_line_x[0] + _B[2] * _delay_line_x[1]
                         + _A[0] * _delay_line_y[0] + _A[1] * _delay_line_y[1];

        _delay_line_x[1] = _delay_line_x[0];
        _delay_line_x[0] = input;
        _delay_line_y[1] = _delay_line_y[0];
        _delay_line_y[0] = output;

        return output;
    }

    void reset() {
        _delay_line_x[0] = _delay_line_x[1] = 0;
        _delay_line_y[0] = _delay_line_y[1] = 0;
    }

private:
    T _B[3]{}; // b0, b1, b2
    T _A[2]{}; // a1, a2

    T _delay_line_x[2]{};   // x[n-1], x[n-2]
    T _delay_line_y[2]{};   // y[n-1], y[n-2]
};

/**
 * @}
 */

}