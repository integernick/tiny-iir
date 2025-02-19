#pragma once

#include "EllipticSolver.h"
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

    void init(double pass_ripple_db, double stop_ripple_db);

    [[nodiscard]] PoleZeroPair get_pole_zero_pairs_s_plane(unsigned int i) final;

    [[nodiscard]] PoleZeroPair get_pole_zero_real_axis() final;

    double _k;
    EllipticSolver _elliptic_solver_k;

    double _k1;
    EllipticSolver _elliptic_solver_k1;

    double _v0; // Solution of sn(1i*v0*N*K1, k1) = j/eps_p
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
    init(pass_ripple_db, stop_ripple_db);
    IIRFilter<N, T, PASS_TYPE>::calculate_cascades(normalized_cutoff_frequency);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
template<FilterPassType _PT, typename>
void
IIRElliptic<N, T, PASS_TYPE>::configure(double normalized_lowcut_freq, double normalized_highcut_freq,
                                        double pass_ripple_db, double stop_ripple_db) {
    init(pass_ripple_db, stop_ripple_db);
    IIRFilter<N, T, PASS_TYPE>::calculate_cascades(normalized_lowcut_freq, normalized_highcut_freq);
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
PoleZeroPair IIRElliptic<N, T, PASS_TYPE>::get_pole_zero_pairs_s_plane(unsigned int i) {
    const double u_i = (2.0 * i + 1.0) / N;
    const double zeta = _elliptic_solver_k.cd(u_i);
    const Complex zero = Complex{0, 1.0 / (_k * zeta)};
    const Complex pole = IMAG_UNIT * _elliptic_solver_k.cd(u_i - IMAG_UNIT * _v0);
    return {pole, zero};
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
PoleZeroPair IIRElliptic<N, T, PASS_TYPE>::get_pole_zero_real_axis() {
    const Complex pole = IMAG_UNIT * _elliptic_solver_k.asn(IMAG_UNIT / _v0);
    const Complex zero = Complex{IIRFilter<N, T, PASS_TYPE>::INFINITY_VALUE, 0};
    return {pole, zero};
}

template<size_t N, typename T, FilterPassType PASS_TYPE>
void IIRElliptic<N, T, PASS_TYPE>::init(double pass_ripple_db, double stop_ripple_db) {
    pass_ripple_db = std::abs(pass_ripple_db);
    stop_ripple_db = std::abs(stop_ripple_db);

    if constexpr (N & 1) {
        IIRFilter<N, T, PASS_TYPE>::_gain_double = 1.0;
    } else {
        IIRFilter<N, T, PASS_TYPE>::_gain_double = std::exp(-pass_ripple_db / 20 * M_LN10);
    }

    const double eps_p = std::sqrt(std::exp(pass_ripple_db * 0.1 * M_LN10) - 1.0);
    const double eps_s = std::sqrt(std::exp(stop_ripple_db * 0.1 * M_LN10) - 1.0);

    const double k1 = eps_p / eps_s;
    _elliptic_solver_k1.init(k1);
    _k = _elliptic_solver_k1.solve_degree_equation(N);
    if (_k > 1.0) {
        return;
    }
    _elliptic_solver_k.init(_k);

    // Solution of sn(j*v0*N*K1, k1) = j/eps_p
    Complex v0 = -IMAG_UNIT * _elliptic_solver_k1.asn(IMAG_UNIT / eps_p) / static_cast<double>(N);
    _v0 = v0.real();
}

}