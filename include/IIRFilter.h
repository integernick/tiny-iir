#pragma once

#include <CascadeFilter.h>

#include <complex>

namespace tiny_iir {

enum class FilterPassType {
    LOW_PASS,
    HIGH_PASS,
};

struct PoleZeroPair {
    std::complex<double> pole;
    std::complex<double> zero;
};

template<size_t N, typename T = double, FilterPassType PASS_TYPE = FilterPassType::LOW_PASS>
class IIRFilter {
public:
    T process(T x) {
        return _cascade_filter.process(x);
    }

    T process(T *x, size_t num_samples) {
        return _cascade_filter.process(x, num_samples);
    }

    [[nodiscard]] const CascadeFilter<N, T> &cascade_filter() const {
        return _cascade_filter;
    }

    virtual PoleZeroPair get_pole_zero_pairs_s_plane(unsigned int i) = 0;

    virtual PoleZeroPair get_pole_zero_real_axis() = 0;

protected:
    virtual ~IIRFilter() = default;

    void configure(double normalized_cutoff_frequency) {
        _cascade_filter.reset();

        if constexpr (PASS_TYPE == FilterPassType::LOW_PASS) {
            _wn_warped = std::tan((M_PI * normalized_cutoff_frequency) / 2);
        } else if constexpr (PASS_TYPE == FilterPassType::HIGH_PASS) {
            _wn_warped = 1.0 / std::tan((M_PI * normalized_cutoff_frequency) / 2);
        }

        // Start from the lowest Q factor pole (closest to the real axis)
        if constexpr (N & 1) {
            add_pole_zero_single_pair(get_pole_zero_real_axis());
        }
        for (int i = static_cast<int>(N) / 2 - 1; i >= 0; --i) {
            add_pole_zero_conjugates_pair(get_pole_zero_pairs_s_plane(i));
        }

        _cascade_filter.init_biquad_cascades();
    }

    static constexpr double INFINITY_VALUE = std::numeric_limits<double>::infinity();

    CascadeFilter<N, T> _cascade_filter;

private:
    std::complex<double> transform(std::complex<double> s) {
        if constexpr (PASS_TYPE == FilterPassType::LOW_PASS) {
            if (s.real() != INFINITY_VALUE) {
                return (1.0 + s) / (1.0 - s);
            } else {
                return {-1, 0};
            }
        } else if constexpr (PASS_TYPE == FilterPassType::HIGH_PASS) {
            if (s.real() != INFINITY_VALUE) {
                return -(1.0 + s) / (1.0 - s);
            } else {
                return {1, 0};
            }
        }
    }

    void add_biquad_coefficients(double b0, double b1, double b2, double a1, double a2) {
        double biquad_gain = 1.0;
        if constexpr (PASS_TYPE == FilterPassType::LOW_PASS) {
            // Evaluate at z = 1 (s = 0)
            biquad_gain = (1 + b1 + b2) / (1 + a1 + a2);
        } else if constexpr (PASS_TYPE == FilterPassType::HIGH_PASS) {
            // Evaluate at z = -1 (s = infinity)
            biquad_gain = (1 - b1 + b2) / (1 - a1 + a2);
        }
        _cascade_filter.push_biquad_coefficients(b0, b1, b2, a1, a2);
        _cascade_filter.set_gain(_cascade_filter.get_gain() / biquad_gain);
    }

    void add_pole_zero_conjugates_pair(const PoleZeroPair &pole_zero_pair) {
        const std::complex<double> pole_s = _wn_warped * pole_zero_pair.pole;
        const std::complex<double> zero_s = _wn_warped * pole_zero_pair.zero;

        const std::complex<double> pole_z = transform(pole_s);
        const std::complex<double> zero_z = transform(zero_s);

        // (z - z0)(z - z0*)=z^2 - z(z0+z0*) + |z0|^2 = z^2 * (1 - 2*Re(z0) * z^-1 + |z0|^2 * z^-2)
        double b0 = 1.0;
        double b1 = -2.0 * zero_z.real();
        double b2 = zero_z.real() * zero_z.real() + zero_z.imag() * zero_z.imag();
        double a1 = -2.0 * pole_z.real();
        double a2 = pole_z.real() * pole_z.real() + pole_z.imag() * pole_z.imag();

        add_biquad_coefficients(b0, b1, b2, a1, a2);
    }

    void add_pole_zero_single_pair(const PoleZeroPair &pole_zero_pair) {
        const std::complex<double> pole_s = _wn_warped * pole_zero_pair.pole;
        const std::complex<double> zero_s = _wn_warped * pole_zero_pair.zero;

        const std::complex<double> pole_z = transform(pole_s);
        const std::complex<double> zero_z = transform(zero_s);

        // (z - z0) = z * (1 - z^-1 * z0 + 0)
        double b0 = 1.0;
        double b1 = -zero_z.real();
        double b2 = 0;
        double a1 = -pole_z.real();
        double a2 = 0;

        add_biquad_coefficients(b0, b1, b2, a1, a2);
    }

    double _wn_warped = 0;  // Wn_warped over 2, technically
};

}