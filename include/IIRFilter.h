#pragma once

#include <CascadeFilter.h>
#include <cmath>
#include <iostream>
#include <complex>

namespace fast_iir {

struct PoleZeroConjugatePair {
    std::complex<double> pole;
    std::complex<double> zero;
};

template<size_t N, typename T = double>
class IIRFilterLPF {
public:
    virtual ~IIRFilterLPF() = default;

    [[nodiscard]] const CascadeFilter<N, T> &cascade_filter() const {
        return _cascade_filter;
    }

    virtual PoleZeroConjugatePair get_pole_zero_pairs_s_plane(unsigned int i) = 0;

protected:
    std::complex<double> transform(std::complex<double> s) {
        if (s.real() != INFINITY_VALUE) {
            return -(s + 1.0) / (s - 1.0);
        } else {
            return {-1, 0};
        }
    }

    void configure_poles_zeros(double normalized_cutoff_frequency) {
        _wn_warped = std::tan((M_PI * normalized_cutoff_frequency) / 2);

        for (int i = 0; i < (N + 1) / 2; ++i) {
            PoleZeroConjugatePair pole_zero_pair = get_pole_zero_pairs_s_plane(i);
            std::cout << "i: " << i << std::endl;
            std::cout << "pole: " << pole_zero_pair.pole << std::endl;
            std::cout << "zero: " << pole_zero_pair.zero << std::endl;

            const std::complex<double> &pole_s = _wn_warped * pole_zero_pair.pole;
            const std::complex<double> &zero_s = _wn_warped * pole_zero_pair.zero;

            const std::complex<double> pole_z = transform(pole_s);
            const std::complex<double> zero_z = transform(zero_s);

            double b0, b1, b2, a1, a2;

            constexpr bool ODD_ORDER = (N & 1);
            if constexpr (ODD_ORDER) {
                if (i != N / 2) {
                    _gain *= (1.0 - pole_z);
                    _gain *= (1.0 - std::conj(pole_z));

                    b0 = 1;
                    b1 = -2 * zero_z.real();
                    b2 = zero_z.real() * zero_z.real() + zero_z.imag() * zero_z.imag();

                    a1 = -2 * pole_z.real();
                    a2 = pole_z.real() * pole_z.real() + pole_z.imag() * pole_z.imag();
                } else {
                    // Single pole on the real axis
                    _gain *= (1.0 - pole_z.real());
                    b0 = 1;
                    b1 = -zero_z.real();
                    b2 = 0;
                    a1 = -pole_z.real();
                    a2 = 0;
                }
            } else {
                _gain *= (1.0 - pole_z);
                _gain *= (1.0 - std::conj(pole_z));
                a1 = -2 * pole_z.real();
                a2 = pole_z.real() * pole_z.real() + pole_z.imag() * pole_z.imag();
            }

            _cascade_filter.push_biquad_coefficients(b0, b1, b2, a1, a2);
        }

        const double k = _gain.real() / (1 << N);
        _cascade_filter.set_gain(k);
    }

    static constexpr double INFINITY_VALUE = std::numeric_limits<double>::infinity();

    CascadeFilter<N, T> _cascade_filter;
    double _wn_warped = 0;  // Wn_warped over 2, technically
    std::complex<double> _gain = 1.0;
};

}