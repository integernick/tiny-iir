#pragma once

#include "utils.h"
#include <CascadeFilter.h>

namespace tiny_iir {

template<size_t N, typename T = double, FilterPassType PASS_TYPE = FilterPassType::LOW_PASS>
class IIRFilter {
public:
    using ValueType = T;

    template<FilterPassType>
    class PassTypeData;

    template<>
    class PassTypeData<FilterPassType::LOW_PASS> {
    public:
        void init(double normalized_cutoff_frequency) {
            normalized_cutoff_frequency = constrain(std::abs(normalized_cutoff_frequency), 0.0, 1.0);
            _wn_warped = std::tan((M_PI * normalized_cutoff_frequency) * 0.5);
        }

        Complex transform(Complex s) {
            if (s.real() == INFINITY_VALUE) {
                return {-1, 0};
            }

            s *= _wn_warped;
            return (1.0 + s) / (1.0 - s);
        }

        [[nodiscard]] double calculate_gain(double b0, double b1, double b2, double a1, double a2) {
            // Evaluate at z = 1 (s = 0)
            return (b0 + b1 + b2) / (1 + a1 + a2);
        }

        static constexpr size_t CASCADE_ORDER = N;

    private:
        double _wn_warped;
    };

    template<>
    class PassTypeData<FilterPassType::HIGH_PASS> {
    public:
        void init(double normalized_cutoff_frequency) {
            normalized_cutoff_frequency = constrain(std::abs(normalized_cutoff_frequency), 0.0, 1.0);
            _wn_warped = 1.0 / std::tan((M_PI * normalized_cutoff_frequency) * 0.5);
        }

        Complex transform(Complex s) {
            if (s.real() == INFINITY_VALUE) {
                return {1, 0};
            }

            s *= _wn_warped;
            return -(1.0 + s) / (1.0 - s);
        }

        [[nodiscard]] double calculate_gain(double b0, double b1, double b2, double a1, double a2) {
            // Evaluate at z = -1 (s = infinity)
            return (b0 - b1 + b2) / (1 - a1 + a2);
        }

        static constexpr size_t CASCADE_ORDER = N;

    private:
        double _wn_warped;
    };

    template<>
    struct PassTypeData<FilterPassType::BAND_PASS> {
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

        [[nodiscard]] double calculate_gain(double b0, double b1, double b2, double a1, double a2) {
            // Evaluate at z = e^(j*wn) (at the peak frequency)
            const Complex z_inv = std::polar(1.0, -_wn_peak);
            const Complex z_inv2 = std::polar(1.0, -2 * _wn_peak);
            Complex response = b0 + b1 * z_inv + b2 * z_inv2;
            response /= (1.0 + a1 * z_inv + a2 * z_inv2);

            return std::abs(response);
        }

        static constexpr size_t CASCADE_ORDER = 2 * N;

    private:
        double _l_inv;
        double _wn_peak;
        double _alpha;
        double _alpha2;
    };

    template<>
    struct PassTypeData<FilterPassType::BAND_STOP> {
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
            _mapped_infinity_value = {(1 - w0_sq) / (1 + w0_sq),
                                      2 * _w0 / (1 + w0_sq)};
        }

        std::pair<Complex, Complex> transform(Complex s) {
            if (s.real() == INFINITY_VALUE) {
                return {_mapped_infinity_value, std::conj(_mapped_infinity_value)};
            } else {
                const double w0_sq = _w0 * _w0;
                const double bw_sq = _bw * _bw;
                const Complex sqrtD = std::sqrt(bw_sq - 4.0 * w0_sq * s * s);
                const Complex den = (s * (1 + w0_sq) - _bw);
                const Complex B = s * (1 - w0_sq);
                const Complex z1 = (B + sqrtD) / den;
                const Complex z2 = (B - sqrtD) / den;
                return {z1, z2};
            }
        }

        [[nodiscard]] double calculate_gain(double b0, double b1, double b2, double a1, double a2) {
            // Evaluate at z = 1 (s = 0)
            return (b0 + b1 + b2) / (1 + a1 + a2);
        }

        static constexpr size_t CASCADE_ORDER = 2 * N;

    private:
        double _w0, _bw;
        Complex _mapped_infinity_value;
    };


    template<typename U>
    [[nodiscard]] T process(U x);

    void process(const T *x, T *out, size_t num_samples);

    [[nodiscard]] T process(const T *x, size_t num_samples);

    template<typename U, std::enable_if_t<std::is_same_v<U, double>, int> = 0>
    void process(const U *x, T *out, size_t num_samples);

    template<typename U, std::enable_if_t<std::is_same_v<U, double>, int> = 0>
    [[nodiscard]] T process(const U *x, size_t num_samples);

    [[nodiscard]] T get_gain() const;

    [[nodiscard]] const T *get_coefficients() const;

    [[nodiscard]] uint32_t get_number_of_blocks() const;

    void print_coefficients() const;

    void reset();

    virtual PoleZeroPair get_pole_zero_pairs_s_plane(unsigned int i) = 0;

    virtual PoleZeroPair get_pole_zero_real_axis() = 0;

    static constexpr size_t ORDER = N;

protected:
    virtual ~IIRFilter() = default;

    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::LOW_PASS
                                         || _PT == FilterPassType::HIGH_PASS)>>
    void calculate_cascades(double normalized_cutoff_frequency) {
        _cascade_filter.reset();

        _pass_type_data.init(normalized_cutoff_frequency);

        // Start from the lowest Q factor pole (closest to the real axis)
        if constexpr (N & 1) {
            PoleZeroPair pole_zero_pair = get_pole_zero_real_axis();
            pole_zero_pair.pole = _pass_type_data.transform(pole_zero_pair.pole);
            pole_zero_pair.zero = _pass_type_data.transform(pole_zero_pair.zero);
            add_pole_zero_single_pair(pole_zero_pair);
        }
        for (int i = static_cast<int>(N) / 2 - 1; i >= 0; --i) {
            PoleZeroPair pole_zero_pair = get_pole_zero_pairs_s_plane(i);
            pole_zero_pair.pole = _pass_type_data.transform(pole_zero_pair.pole);
            pole_zero_pair.zero = _pass_type_data.transform(pole_zero_pair.zero);
            add_pole_zero_conjugates_pair(pole_zero_pair);
        }

        _cascade_filter.set_gain(_gain_double);
        _cascade_filter.init_biquad_cascades();
    }

    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::BAND_PASS
                                         || _PT == FilterPassType::BAND_STOP)>>
    void calculate_cascades(double cutlow_freq, double cuthigh_freq) {
        _cascade_filter.reset();

        constexpr double MIN_FREQ = 1e-8;
        constexpr double MAX_FREQ = 2 * M_PI - MIN_FREQ;
        if (cutlow_freq < MIN_FREQ) {
            cutlow_freq = MIN_FREQ;
        }
        if (cuthigh_freq > MAX_FREQ) {
            cuthigh_freq = MAX_FREQ;
        }

        _pass_type_data.init(cutlow_freq, cuthigh_freq);

        auto add_two_pole_zero_pairs_from_one = [&](const PoleZeroPair &pole_zero_pair) {
            auto [pole1, pole2] = _pass_type_data.transform(pole_zero_pair.pole);
            auto [zero1, zero2] = _pass_type_data.transform(pole_zero_pair.zero);

            PoleZeroPair pz1(pole1, zero1);
            add_pole_zero_conjugates_pair(pz1);

            PoleZeroPair pz2(pole2, zero2);
            add_pole_zero_conjugates_pair(pz2);
        };

        if constexpr (N & 1) {
            const PoleZeroPair pole_zero_pair = get_pole_zero_real_axis();
            auto [pole1, pole2] = _pass_type_data.transform(pole_zero_pair.pole);
            auto [zero1, zero2] = _pass_type_data.transform(pole_zero_pair.zero);

            PoleZeroPair pz1(pole1, zero1);
            PoleZeroPair pz2(pole2, zero2);
            add_pole_zero_pairs({pz1, pz2});
        }
        // Start from the lowest Q factor pole (closest to the real axis)
        for (int i = static_cast<int>(N) / 2 - 1; i >= 0; --i) {
            PoleZeroPair pole_zero_pair = get_pole_zero_pairs_s_plane(i);
            add_two_pole_zero_pairs_from_one(pole_zero_pair);
        }

        _cascade_filter.set_gain(_gain_double);
        _cascade_filter.init_biquad_cascades();
    }

    static constexpr double INFINITY_VALUE = std::numeric_limits<double>::infinity();

    CascadeFilter<PassTypeData<PASS_TYPE>::CASCADE_ORDER, T> _cascade_filter;
    double _gain_double = 1.0;

private:
    void add_biquad_coefficients(double b0, double b1, double b2, double a1, double a2);

    void add_pole_zero_conjugates_pair(const PoleZeroPair &pole_zero_pair);

    void add_pole_zero_single_pair(const PoleZeroPair &pole_zero_pair);

    void add_pole_zero_pairs(const std::pair<PoleZeroPair, PoleZeroPair> &pole_zero_pairs);

    PassTypeData<PASS_TYPE> _pass_type_data;
};

template<size_t N, typename T, FilterPassType PASS_TYPE>
T IIRFilter<N, T, PASS_TYPE>::get_gain() const {
    return _cascade_filter.get_gain();
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
const T *IIRFilter<N, T, PASS_TYPE>::get_coefficients() const {
    return _cascade_filter.get_coefficients();
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
uint32_t IIRFilter<N, T, PASS_TYPE>::get_number_of_blocks() const {
    return _cascade_filter.get_number_of_blocks();
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
void IIRFilter<N, T, PASS_TYPE>::print_coefficients() const {
    _cascade_filter.print_coefficients();
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
T IIRFilter<N, T, PASS_TYPE>::process(const T *x, size_t num_samples) {
    return _cascade_filter.process(x, num_samples);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
template<typename U, std::enable_if_t<std::is_same_v<U, double>, int>>
void IIRFilter<N, T, PASS_TYPE>::process(const U *x, T *out, size_t num_samples) {
    _cascade_filter.process(x, out, num_samples);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
template<typename U, std::enable_if_t<std::is_same_v<U, double>, int>>
T IIRFilter<N, T, PASS_TYPE>::process(const U *x, size_t num_samples) {
    return _cascade_filter.process(x, num_samples);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
void IIRFilter<N, T, PASS_TYPE>::reset() {
    _cascade_filter.reset_state();
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
void IIRFilter<N, T, PASS_TYPE>::add_biquad_coefficients(double b0, double b1, double b2, double a1, double a2) {
    _cascade_filter.push_biquad_coefficients(b0, b1, b2, a1, a2);
    const double biquad_gain = _pass_type_data.calculate_gain(b0, b1, b2, a1, a2);
    _gain_double /= biquad_gain;
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
void IIRFilter<N, T, PASS_TYPE>::add_pole_zero_conjugates_pair(const PoleZeroPair &pole_zero_pair) {
    const Complex &pole_z = pole_zero_pair.pole;
    const Complex &zero_z = pole_zero_pair.zero;

    // (z - z0)(z - z0*)=z^2 - z(z0+z0*) + |z0|^2 = z^2 * (1 - 2*Re(z0) * z^-1 + |z0|^2 * z^-2)
    double b0 = 1.0;
    double b1 = -2.0 * zero_z.real();
    double b2 = zero_z.real() * zero_z.real() + zero_z.imag() * zero_z.imag();
    double a1 = -2.0 * pole_z.real();
    double a2 = pole_z.real() * pole_z.real() + pole_z.imag() * pole_z.imag();

    add_biquad_coefficients(b0, b1, b2, a1, a2);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
void IIRFilter<N, T, PASS_TYPE>::add_pole_zero_single_pair(const PoleZeroPair &pole_zero_pair) {
    const Complex &pole_z = pole_zero_pair.pole;
    const Complex &zero_z = pole_zero_pair.zero;

    // (z - z0) = z * (1 - z^-1 * z0 + 0)
    double b0 = 1.0;
    double b1 = -zero_z.real();
    double b2 = 0;
    double a1 = -pole_z.real();
    double a2 = 0;

    add_biquad_coefficients(b0, b1, b2, a1, a2);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
void IIRFilter<N, T, PASS_TYPE>::add_pole_zero_pairs(const std::pair<PoleZeroPair, PoleZeroPair> &pole_zero_pairs) {
    // (z - z1)(z - z2)=z^2 - z(z1+z2) + z1z2 = z^2 * (1 - (z1+z2)z^-1 + (z1*z2)*z^-2)
    const Complex &p1 = pole_zero_pairs.first.pole;
    const Complex &p2 = pole_zero_pairs.second.pole;
    const Complex &z1 = pole_zero_pairs.first.zero;
    const Complex &z2 = pole_zero_pairs.second.zero;

    double b0 = 1.0;
    double b1 = -(pole_zero_pairs.first.zero.real() + pole_zero_pairs.second.zero.real());
    double b2 = z1.imag() == 0
                ? z1.real() * z2.real()
                : std::abs(z1);

    double a1 = -(pole_zero_pairs.first.pole.real() + pole_zero_pairs.second.pole.real());
    double a2 = p1.imag() == 0
                ? p1.real() * p2.real()
                : std::norm(p1);

    add_biquad_coefficients(b0, b1, b2, a1, a2);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
template<typename U>
T IIRFilter<N, T, PASS_TYPE>::process(U x) {
    return _cascade_filter.process(x);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
void IIRFilter<N, T, PASS_TYPE>::process(const T *x, T *out, size_t num_samples) {
    _cascade_filter.process(x, out, num_samples);
}

}