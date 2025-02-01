#include "test_utils.h"
#include "test_utils_fixed_point.h"

#include <butter/IIRButter.h>

using namespace tiny_iir;

/* TODO: Check why MATLAB calculates b0, b1, b2 coefficients that slightly differ from 1.0 or 2.0 by about 0.01
   Therefore the impulse response is also slightly different! */
constexpr double TOL = 1e-4;

TEST(ButterTest, ButterLPFQ31Coeffs) {
    IIRButter<7, q31_t> butter_lpf(0.75);

    butter_lpf.cascade_filter().print_coefficients();

    constexpr double GAIN_EXPECTED = 0.161071111239757;
    // TODO: Check why MATLAB calculates b0, b1, b2 coefficients that slightly differ from 1.0 or 2.0 by about 0.01
    std::vector<double> expected_coeffs = {
            1, 1.0, 0.0, -0.414213562373084, -0.0,
            1, 2.0, 1.0, -0.863862810051667, -0.221686502004782,
            1, 2.0, 1.0, -0.981497128142559, -0.388046550049447,
            1, 2.0, 1.0, -1.22194528464162, -0.728091594018049,
    };
    normalize_coeffs(expected_coeffs);

    test_coeffs(butter_lpf, GAIN_EXPECTED, expected_coeffs, TOL);
}

TEST(ButterTest, ButterLPFQ31ImpulseResponse) {
    IIRButter<12, q31_t> butter_lpf(0.3);

    const std::vector<double> expected = {
            6.86125789668954e-06, 0.000115198364232461, 0.000924704457530324, 0.00472285891438076,
            0.0172210594453799, 0.0476308528941755, 0.103457117533882, 0.17971446218948,
            0.250511907492351, 0.275449809620748, 0.224243365851443, 0.104332560457353,
            -0.033361970739935, -0.121822893976523, -0.121804233594017, -0.0487733669399159,
            0.0390135578632716, 0.0824690366122271, 0.0613263110554832, 0.00248849970475534
    };

    test_impulse_response(butter_lpf, expected, TOL);
    butter_lpf.reset();
    test_impulse_response_batch(butter_lpf, expected, TOL);
}
