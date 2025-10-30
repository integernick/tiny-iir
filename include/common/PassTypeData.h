#pragma once

#include "common_utils.h"

#include <numbers>

namespace tiny_iir {
template<FilterPassType, uint32_t N, typename DESIGN_T>
class PassTypeData;

template<uint32_t N, typename DESIGN_T>
class PassTypeData<FilterPassType::LowPass, N, DESIGN_T> {
    static_assert(std::is_same_v<DESIGN_T, double> or std::is_same_v<DESIGN_T, float>,
                  "DESIGN_T must be float or double");

public:
    void init(DESIGN_T normalized_cutoff_frequency) {
        normalized_cutoff_frequency
                = constrain(std::abs(normalized_cutoff_frequency), DESIGN_T{0}, DESIGN_T{1});
        _wn_warped = std::tan(
                (std::numbers::pi_v<DESIGN_T> * normalized_cutoff_frequency) * static_cast<DESIGN_T>(0.5));
    }

    Complex<DESIGN_T> transform(Complex<DESIGN_T> s) {
        if (s.real() == std::numeric_limits<DESIGN_T>::infinity()) {
            return {-DESIGN_T{1}, 0};
        }

        s *= _wn_warped;
        return (DESIGN_T{1} + s) / (DESIGN_T{1} - s);
    }

    [[nodiscard]] double calculate_gain(const BiquadCoefficients<DESIGN_T> &biquad_coefficients) {
        // Evaluate at z = 1 (s = 0)
        return (biquad_coefficients.b0 + biquad_coefficients.b1 + biquad_coefficients.b2)
               / (DESIGN_T{1} + biquad_coefficients.a1 + biquad_coefficients.a2);
    }

    static constexpr uint32_t CASCADE_ORDER = N;

private:
    double _wn_warped;
};

template<uint32_t N, typename DESIGN_T>
class PassTypeData<FilterPassType::HighPass, N, DESIGN_T> {
    static_assert(std::is_same_v<DESIGN_T, double> or std::is_same_v<DESIGN_T, float>,
                  "DESIGN_T must be float or double");

public:
    void init(DESIGN_T normalized_cutoff_frequency) {
        normalized_cutoff_frequency
                = constrain(std::abs(normalized_cutoff_frequency), DESIGN_T{0}, DESIGN_T{1});
        _wn_warped = DESIGN_T{1} / std::tan(
                (std::numbers::pi_v<DESIGN_T> * normalized_cutoff_frequency) * static_cast<DESIGN_T>(0.5));
    }

    Complex<DESIGN_T> transform(Complex<DESIGN_T> s) {
        if (s.real() == std::numeric_limits<DESIGN_T>::infinity()) {
            return {1, 0};
        }

        s *= _wn_warped;
        return -(DESIGN_T{1} + s) / (DESIGN_T{1} - s);
    }

    [[nodiscard]] DESIGN_T calculate_gain(const BiquadCoefficients<DESIGN_T> &biquad_coefficients) {
        // Evaluate at z = -1 (s = infinity)
        return (biquad_coefficients.b0 - biquad_coefficients.b1 + biquad_coefficients.b2)
               / (DESIGN_T{1} - biquad_coefficients.a1 + biquad_coefficients.a2);
    }

    static constexpr uint32_t CASCADE_ORDER = N;

private:
    DESIGN_T _wn_warped;
};

template<uint32_t N, typename DESIGN_T>
class PassTypeData<FilterPassType::BandPass, N, DESIGN_T> {
    static_assert(std::is_same_v<DESIGN_T, double> or std::is_same_v<DESIGN_T, float>,
                  "DESIGN_T must be float or double");

public:
    void init(DESIGN_T cutlow_freq, DESIGN_T cuthigh_freq) {
        static constexpr DESIGN_T m_pi_2_v = std::numbers::pi_v<DESIGN_T> * static_cast<DESIGN_T>(0.5);

        cutlow_freq = constrain(std::abs(cutlow_freq), 0.0, DESIGN_T{1});
        cuthigh_freq = constrain(std::abs(cuthigh_freq), 0.0, DESIGN_T{1});

        if (cutlow_freq > cuthigh_freq) {
            DESIGN_T temp = cuthigh_freq;
            cuthigh_freq = cutlow_freq;
            cutlow_freq = temp;
        }

        // s = l * (z^2 - 2*a*z + 1) / (z^2 - 1)
        const DESIGN_T w_central = m_pi_2_v * (cuthigh_freq + cutlow_freq) * 0.5;
        const DESIGN_T w_bandwidth = m_pi_2_v * (cuthigh_freq - cutlow_freq);
        const DESIGN_T w1 = w_central - w_bandwidth * 0.5;
        const DESIGN_T w2 = w_central + w_bandwidth * 0.5;

        _l_inv = std::tan(m_pi_2_v * (cuthigh_freq - cutlow_freq));
        _alpha = std::cos(m_pi_2_v * (cuthigh_freq + cutlow_freq))
                 / std::cos(m_pi_2_v * (cuthigh_freq - cutlow_freq));
        _alpha2 = _alpha * _alpha;
        const DESIGN_T w_avg = std::sqrt(std::tan(w1) * std::tan(w2));
        _wn_peak = 2 * std::atan(w_avg);
    }

    std::pair<Complex<DESIGN_T>, Complex<DESIGN_T>> transform(Complex<DESIGN_T> s) {
        if (s.real() == std::numeric_limits<DESIGN_T>::infinity()) {
            return {{-DESIGN_T{1}, 0},
                    {DESIGN_T{1},  0}};
        }

        s *= _l_inv;
        const Complex<DESIGN_T> A = DESIGN_T{1} - s;
        const Complex<DESIGN_T> C = DESIGN_T{1} + s;
        const Complex<DESIGN_T> sqrtD = std::sqrt(_alpha2 - A * C);
        std::complex z1 = (_alpha + sqrtD) / A;
        std::complex z2 = (_alpha - sqrtD) / A;

        return {z1, z2};
    }

    [[nodiscard]] DESIGN_T
    calculate_gain(const BiquadCoefficients<DESIGN_T> &biquad_coefficients) {
        // Evaluate at z = e^(j*wn) (at the peak frequency)
        const Complex<DESIGN_T> z_inv = std::polar(DESIGN_T{1}, -_wn_peak);
        const Complex<DESIGN_T> z_inv2 = std::polar(DESIGN_T{1}, -2 * _wn_peak);
        Complex<DESIGN_T> response = biquad_coefficients.b0 + biquad_coefficients.b1 * z_inv + biquad_coefficients.b2 * z_inv2;
        response /= (DESIGN_T{1} + biquad_coefficients.a1 * z_inv + biquad_coefficients.a2 * z_inv2);

        return std::abs(response);
    }

    static constexpr uint32_t CASCADE_ORDER = 2 * N;

private:
    DESIGN_T _l_inv;
    DESIGN_T _wn_peak;
    DESIGN_T _alpha;
    DESIGN_T _alpha2;
};

template<uint32_t N, typename DESIGN_T>
class PassTypeData<FilterPassType::BandStop, N, DESIGN_T> {
    static_assert(std::is_same_v<DESIGN_T, double> or std::is_same_v<DESIGN_T, float>,
                  "DESIGN_T must be float or double");

public:
    void init(DESIGN_T cutlow_freq, DESIGN_T cuthigh_freq) {
        static constexpr DESIGN_T m_pi_2_v = std::numbers::pi_v<DESIGN_T> * static_cast<DESIGN_T>(0.5);

        cutlow_freq = constrain(std::abs(cutlow_freq), 0.0, DESIGN_T{1});
        cuthigh_freq = constrain(std::abs(cuthigh_freq), 0.0, DESIGN_T{1});

        if (cutlow_freq > cuthigh_freq) {
            DESIGN_T temp = cuthigh_freq;
            cuthigh_freq = cutlow_freq;
            cutlow_freq = temp;
        }

        DESIGN_T w1n = std::tan(m_pi_2_v * cutlow_freq);
        DESIGN_T w2n = std::tan(m_pi_2_v * cuthigh_freq);
        _bw = w2n - w1n;
        _w0 = std::sqrt(w1n * w2n);

        const DESIGN_T w0_sq = _w0 * _w0;
        _mapped_infinity_value
                = {(1 - w0_sq) / (1 + w0_sq), 2 * _w0 / (1 + w0_sq)};
    }

    std::pair<Complex<DESIGN_T>, Complex<DESIGN_T>> transform(Complex<DESIGN_T> s) {
        if (s.real() == std::numeric_limits<DESIGN_T>::infinity()) {
            return {_mapped_infinity_value, std::conj(_mapped_infinity_value)};
        }

        const DESIGN_T w0_sq = _w0 * _w0;
        const DESIGN_T bw_sq = _bw * _bw;
        const Complex<DESIGN_T> sqrtD = std::sqrt(bw_sq - 4.0 * w0_sq * s * s);
        const Complex<DESIGN_T> den = (s * (1 + w0_sq) - _bw);
        const Complex<DESIGN_T> B = s * (1 - w0_sq);
        const Complex<DESIGN_T> z1 = (B + sqrtD) / den;
        const Complex<DESIGN_T> z2 = (B - sqrtD) / den;
        return {z1, z2};
    }

    [[nodiscard]] DESIGN_T calculate_gain(const BiquadCoefficients<DESIGN_T> &biquad_coefficients) {
        // Evaluate at z = 1 (s = 0)
        return (biquad_coefficients.b0 + biquad_coefficients.b1 + biquad_coefficients.b2)
               / (DESIGN_T{1} + biquad_coefficients.a1 + biquad_coefficients.a2);
    }

    static constexpr uint32_t CASCADE_ORDER = 2 * N;

private:
    DESIGN_T _w0;
    DESIGN_T _bw;
    Complex<DESIGN_T> _mapped_infinity_value;
};
} // namespace tiny_iir