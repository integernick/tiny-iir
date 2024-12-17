#pragma once

#include <IIRFilter.h>

namespace fast_iir {
    template<size_t N = 2, typename T = double>
    class Cheby1LPF : public IIRFilterLPF<N, T> {
    public:
        Cheby1LPF(double normalized_cutoff_frequency, double ripple_db) {
            configure(normalized_cutoff_frequency, ripple_db);
        }

        PoleZeroConjugatePair get_pole_zero_pairs_s_plane(unsigned int i) {
            const double phi = (2 * i + 1) * _d_phi;
            const double sin_phi = std::sin(phi);
            const double cos_phi = std::cos(phi);
            const double pole_s_real = -_sinh_mu * sin_phi;
            const double pole_s_imag = _cosh_mu * cos_phi;
            return {{pole_s_real, pole_s_imag},
                    {IIRFilterLPF<N, T>::INFINITY_VALUE, 0}};
        }

        void configure(double Wn, double ripple_db) {
            if constexpr (N & 1) {
                IIRFilterLPF<N, T>::_gain = 1.0;
            } else {
                IIRFilterLPF<N, T>::_gain = std::exp(-ripple_db / 20 * M_LN10);
            }

            const double epsilon = std::sqrt(std::exp(ripple_db * 0.1 * M_LN10) - 1);
            const double mu = std::asinh(1.0 / epsilon) / N;
            _d_phi = M_PI_2 / N;
            _sinh_mu = std::sinh(mu);
            _cosh_mu = std::cosh(mu);

            IIRFilterLPF<N, T>::configure_poles_zeros(Wn);
        }

    private:
        double _d_phi = 0;
        double _sinh_mu = 0;
        double _cosh_mu = 0;
    };

}