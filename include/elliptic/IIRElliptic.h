#pragma once

#include "elliptic_utils.h"
#include "common/IIRFilter.h"

namespace tiny_iir {

/**
 * @brief Elliptic filter
 *
 * @tparam N   Filter order.
 * @tparam T   Data type.
 * @tparam PASS_TYPE   Pass type (low, high).
 */
template<size_t N = 2, typename T = double, FilterPassType PASS_TYPE = FilterPassType::LOW_PASS>
class IIRElliptic : public IIRFilter<N, T, PASS_TYPE> {
public:
    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::LOW_PASS
                                         || _PT == FilterPassType::HIGH_PASS)>>
    IIRElliptic(double normalized_cutoff_frequency, double pass_ripple_db, double stop_ripple_db);

    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::BAND_PASS
                                         || _PT == FilterPassType::BAND_STOP)>>
    IIRElliptic(double normalized_lowcut_freq, double normalized_highcut_freq,
                double pass_ripple_db, double stop_ripple_db);

    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::LOW_PASS
                                         || _PT == FilterPassType::HIGH_PASS)>>
    void configure(double normalized_cutoff_frequency, double pass_ripple_db, double stop_ripple_db);

    template<FilterPassType _PT = PASS_TYPE,
            typename = std::enable_if_t<(_PT == FilterPassType::BAND_PASS
                                         || _PT == FilterPassType::BAND_STOP)>>
    void configure(double normalized_lowcut_freq, double normalized_highcut_freq,
                   double pass_ripple_db, double stop_ripple_db);

private:
    static constexpr size_t L = N / 2;
    static constexpr Complex IMAG_UNIT = Complex{0, 1};

    void init_analog() final;

    double _pass_ripple_db = 0;
    double _stop_ripple_db = 0;
};


template<size_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType _PT, typename>
IIRElliptic<N, T, PASS_TYPE>::IIRElliptic(double normalized_cutoff_frequency, double pass_ripple_db,
                                          double stop_ripple_db)
        : IIRFilter<N, T, PASS_TYPE>() {
    configure(normalized_cutoff_frequency, pass_ripple_db, stop_ripple_db);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType _PT, typename>
IIRElliptic<N, T, PASS_TYPE>::IIRElliptic(double normalized_lowcut_freq, double normalized_highcut_freq,
                                          double pass_ripple_db, double stop_ripple_db)
        : IIRFilter<N, T, PASS_TYPE>() {
    configure(normalized_lowcut_freq, normalized_highcut_freq, pass_ripple_db, stop_ripple_db);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType _PT, typename>
void IIRElliptic<N, T, PASS_TYPE>::configure(double normalized_cutoff_frequency,
                                             double pass_ripple_db, double stop_ripple_db) {
    pass_ripple_db = std::abs(pass_ripple_db);
    stop_ripple_db = std::abs(stop_ripple_db);
    if (pass_ripple_db != _pass_ripple_db || stop_ripple_db != _stop_ripple_db) {
        _pass_ripple_db = pass_ripple_db;
        _stop_ripple_db = stop_ripple_db;
        init_analog();
    }
    IIRFilter<N, T, PASS_TYPE>::calculate_cascades(normalized_cutoff_frequency);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType _PT, typename>
void
IIRElliptic<N, T, PASS_TYPE>::configure(double normalized_lowcut_freq, double normalized_highcut_freq,
                                        double pass_ripple_db, double stop_ripple_db) {
    pass_ripple_db = std::abs(pass_ripple_db);
    stop_ripple_db = std::abs(stop_ripple_db);
    if (pass_ripple_db != _pass_ripple_db || stop_ripple_db != _stop_ripple_db) {
        _pass_ripple_db = pass_ripple_db;
        _stop_ripple_db = stop_ripple_db;
        init_analog();
    }
    IIRFilter<N, T, PASS_TYPE>::calculate_cascades(normalized_lowcut_freq, normalized_highcut_freq);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
void IIRElliptic<N, T, PASS_TYPE>::init_analog() {
    if constexpr (N & 1) {
        IIRFilter<N, T, PASS_TYPE>::_gain_double = 1.0;
    } else {
        IIRFilter<N, T, PASS_TYPE>::_gain_double = std::exp(-_pass_ripple_db / 20 * M_LN10);
    }

    const double eps_p = std::sqrt(std::exp(_pass_ripple_db * 0.1 * M_LN10) - 1.0);
    const double eps_s = std::sqrt(std::exp(_stop_ripple_db * 0.1 * M_LN10) - 1.0);

    const double k1 = eps_p / eps_s; // The ratio (eps_p / eps_s) << 1

    // Solution of sn(j*v0*N*K1, k1) = j/eps_p
    const double K1 = calculate_elliptic_integral(k1);
    const double K1_prime = calculate_elliptic_integral(get_complimentary(k1));
    const double R1 = K1_prime / K1;
    const Complex v0_c = -IMAG_UNIT * asn(IMAG_UNIT / eps_p, k1, R1) / static_cast<double>(N);
    const double v0 = v0_c.real();

    const double k = solve_degree_equation(N, get_complimentary(k1)); // The ratio (w_p / w_s) < 1
    if (k > 1.0) {
        return;
    }
    const double K = calculate_elliptic_integral(k);
    const double K_prime = calculate_elliptic_integral(get_complimentary(k));
    const double R = K_prime / K;

    if constexpr (N & 1) {
        const Complex pole = IMAG_UNIT * asn(IMAG_UNIT / v0, k, R);
        const Complex zero = Complex{IIRFilter<N, T, PASS_TYPE>::INFINITY_VALUE, 0};
        IIRFilter<N, T, PASS_TYPE>::_analog_pole_zero_pairs[(N + 1) / 2 - 1] = {pole, zero};
    }
    for (int i = 0; i < N / 2; ++i) {
        const double u_i = (2.0 * i + 1.0) / N;
        const double zeta = cd(u_i, k).real();
        const Complex zero = Complex{0, 1.0 / (k * zeta)};
        const Complex pole = IMAG_UNIT * cd(u_i - IMAG_UNIT * v0, k);
        IIRElliptic<N, T, PASS_TYPE>::_analog_pole_zero_pairs[i] = {pole, zero};
    }
}

}