#pragma once

#include "common_utils.h"

namespace tiny_iir {
template<FilterPassType, uint32_t N>
class PassTypeData;

template<uint32_t N>
class PassTypeData<FilterPassType::LOW_PASS, N> {
public:
    void init(double normalized_cutoff_frequency) {
        normalized_cutoff_frequency
                = constrain(std::abs(normalized_cutoff_frequency), 0.0, 1.0);
        _wn_warped = std::tan((M_PI * normalized_cutoff_frequency) * 0.5);
    }

    Complex transform(Complex s) {
        if (s.real() == INFINITY_VALUE) {
            return {-1, 0};
        }

        s *= _wn_warped;
        return (1.0 + s) / (1.0 - s);
    }

    [[nodiscard]] double
    calculate_gain(const BiquadCoefficients &biquad_coefficients) {
        // Evaluate at z = 1 (s = 0)
        return (biquad_coefficients.b0 + biquad_coefficients.b1 + biquad_coefficients.b2)
               / (1.0 + biquad_coefficients.a1 + biquad_coefficients.a2);
    }

    static constexpr uint32_t CASCADE_ORDER = N;

private:
    double _wn_warped;
};

template<uint32_t N>
class PassTypeData<FilterPassType::HIGH_PASS, N> {
public:
    void init(double normalized_cutoff_frequency) {
        normalized_cutoff_frequency
                = constrain(std::abs(normalized_cutoff_frequency), 0.0, 1.0);
        _wn_warped = 1.0 / std::tan((M_PI * normalized_cutoff_frequency) * 0.5);
    }

    Complex transform(Complex s) {
        if (s.real() == INFINITY_VALUE) {
            return {1, 0};
        }

        s *= _wn_warped;
        return -(1.0 + s) / (1.0 - s);
    }

    [[nodiscard]] double
    calculate_gain(const BiquadCoefficients &biquad_coefficients) {
        // Evaluate at z = -1 (s = infinity)
        return (biquad_coefficients.b0 - biquad_coefficients.b1 + biquad_coefficients.b2)
               / (1.0 - biquad_coefficients.a1 + biquad_coefficients.a2);
    }

    static constexpr uint32_t CASCADE_ORDER = N;

private:
    double _wn_warped;
};

template<uint32_t N>
class PassTypeData<FilterPassType::BAND_PASS, N> {
public:
    void init(double cutlow_freq, double cuthigh_freq) {
        cutlow_freq = constrain(std::abs(cutlow_freq), 0.0, 1.0);
        cuthigh_freq = constrain(std::abs(cuthigh_freq), 0.0, 1.0);
        if (cutlow_freq > cuthigh_freq) {
            double temp = cuthigh_freq;
            cuthigh_freq = cutlow_freq;
            cutlow_freq = temp;
        }

        // s = l * (z^2 - 2*a*z + 1) / (z^2 - 1)
        const double w_central = M_PI_2 * (cuthigh_freq + cutlow_freq) * 0.5;
        const double w_bandwidth = M_PI_2 * (cuthigh_freq - cutlow_freq);
        const double w1 = w_central - w_bandwidth * 0.5;
        const double w2 = w_central + w_bandwidth * 0.5;

        _l_inv = std::tan(M_PI_2 * (cuthigh_freq - cutlow_freq));
        _alpha = std::cos(M_PI_2 * (cuthigh_freq + cutlow_freq))
                 / std::cos(M_PI_2 * (cuthigh_freq - cutlow_freq));
        _alpha2 = _alpha * _alpha;
        const double w_avg = std::sqrt(std::tan(w1) * std::tan(w2));
        _wn_peak = 2 * std::atan(w_avg);
    }

    std::pair<Complex, Complex> transform(Complex s) {
        if (s.real() == INFINITY_VALUE) {
            return {{-1, 0},
                    {1,  0}};
        }

        s *= _l_inv;
        const Complex A = 1.0 - s;
        const Complex C = 1.0 + s;
        const Complex sqrtD = std::sqrt(_alpha2 - A * C);
        std::complex z1 = (_alpha + sqrtD) / A;
        std::complex z2 = (_alpha - sqrtD) / A;

        return {z1, z2};
    }

    [[nodiscard]] double
    calculate_gain(const BiquadCoefficients &biquad_coefficients) {
        // Evaluate at z = e^(j*wn) (at the peak frequency)
        const Complex z_inv = std::polar(1.0, -_wn_peak);
        const Complex z_inv2 = std::polar(1.0, -2 * _wn_peak);
        Complex response = biquad_coefficients.b0 + biquad_coefficients.b1 * z_inv + biquad_coefficients.b2 * z_inv2;
        response /= (1.0 + biquad_coefficients.a1 * z_inv + biquad_coefficients.a2 * z_inv2);

        return std::abs(response);
    }

    static constexpr uint32_t CASCADE_ORDER = 2 * N;

private:
    double _l_inv;
    double _wn_peak;
    double _alpha;
    double _alpha2;
};

template<uint32_t N>
class PassTypeData<FilterPassType::BAND_STOP, N> {
public:
    void init(double cutlow_freq, double cuthigh_freq) {
        cutlow_freq = constrain(std::abs(cutlow_freq), 0.0, 1.0);
        cuthigh_freq = constrain(std::abs(cuthigh_freq), 0.0, 1.0);

        if (cutlow_freq > cuthigh_freq) {
            double temp = cuthigh_freq;
            cuthigh_freq = cutlow_freq;
            cutlow_freq = temp;
        }

        double w1n = std::tan(M_PI_2 * cutlow_freq);
        double w2n = std::tan(M_PI_2 * cuthigh_freq);
        _bw = w2n - w1n;
        _w0 = std::sqrt(w1n * w2n);

        const double w0_sq = _w0 * _w0;
        _mapped_infinity_value
                = {(1 - w0_sq) / (1 + w0_sq), 2 * _w0 / (1 + w0_sq)};
    }

    std::pair<Complex, Complex> transform(Complex s) {
        if (s.real() == INFINITY_VALUE) {
            return {_mapped_infinity_value, std::conj(_mapped_infinity_value)};
        }

        const double w0_sq = _w0 * _w0;
        const double bw_sq = _bw * _bw;
        const Complex sqrtD = std::sqrt(bw_sq - 4.0 * w0_sq * s * s);
        const Complex den = (s * (1 + w0_sq) - _bw);
        const Complex B = s * (1 - w0_sq);
        const Complex z1 = (B + sqrtD) / den;
        const Complex z2 = (B - sqrtD) / den;
        return {z1, z2};
    }

    [[nodiscard]] double calculate_gain(const BiquadCoefficients &biquad_coefficients) {
        // Evaluate at z = 1 (s = 0)
        return (biquad_coefficients.b0 + biquad_coefficients.b1 + biquad_coefficients.b2)
               / (1.0 + biquad_coefficients.a1 + biquad_coefficients.a2);
    }

    static constexpr uint32_t CASCADE_ORDER = 2 * N;

private:
    double _w0;
    double _bw;
    Complex _mapped_infinity_value;
};
} // namespace tiny_iir