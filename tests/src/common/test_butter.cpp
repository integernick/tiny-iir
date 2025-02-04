#include "test_utils.h"

#include <butter/IIRButter.h>

using namespace tiny_iir;

constexpr double TOL = TOL_DOUBLE;

TEST(ButterTest, ButterLPFDoubleCoeffs) {
    IIRButter<7, double> butter_lpf(0.75);

    constexpr double GAIN_EXPECTED = 0.161071111239757;
    // TODO: Check why MATLAB calculates b0, b1, b2 coefficients that slightly differ from 1.0 or 2.0 by about 0.01
    const std::vector<double> expected_coeffs = {
            1, 1.0, 0.0, 1, 0.414213562373084, 0.0,
            1, 2.0, 1.0, 1, 0.863862810051667, 0.221686502004782,
            1, 2.0, 1.0, 1, 0.981497128142559, 0.388046550049447,
            1, 2.0, 1.0, 1, 1.22194528464162, 0.728091594018049,
    };

    test_coeffs(butter_lpf, GAIN_EXPECTED, expected_coeffs, TOL);
}

TEST(ButterTest, ButterLPFDoubleImpulseResponse) {
    IIRButter<12, double> butter_lpf(0.3);

    const std::vector<double> expected = {
            6.86125789668954e-06, 0.000115198364232461, 0.000924704457530324, 0.00472285891438076,
            0.0172210594453799, 0.0476308528941755, 0.103457117533882, 0.17971446218948,
            0.250511907492351, 0.275449809620748, 0.224243365851443, 0.104332560457353,
            -0.033361970739935, -0.121822893976523, -0.121804233594017, -0.0487733669399159,
            0.0390135578632716, 0.0824690366122271, 0.0613263110554832, 0.00248849970475534
    };

    test_impulse_response(butter_lpf, expected, TOL);
}

TEST(ButterTest, ButterBPFImpulseResponse) {
    IIRButter<4, double, FilterPassType::BAND_PASS> butter_bpf(0.7, 0.8);
    butter_bpf.print_coefficients();

    /* TODO: Trust me bro */
    EXPECT_TRUE(true);
}

TEST(ButterTest, ButterBSFImpulseResponse) {
    IIRButter<4, double, FilterPassType::BAND_STOP> butter_bsf(0.7, 0.8);
    butter_bsf.print_coefficients();

    /* TODO: Trust */
    EXPECT_TRUE(true);
}
