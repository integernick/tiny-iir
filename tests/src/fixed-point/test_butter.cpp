#include "test_utils.h"
#include "test_utils_fixed_point.h"

#include <butter/IIRButter.h>

using namespace tiny_iir;

/* TODO: Check why MATLAB calculates b0, b1, b2 coefficients that slightly differ from 1.0 or 2.0 by about 0.01
   Therefore the impulse response is also slightly different! */
constexpr double TOL = 1e-4;

TEST(ButterTest, ButterLPFQ31Coeffs) {
    IIRButter<7, q31_t> butter_lpf(0.75);
    constexpr double GAIN_EXPECTED = 0.161071111239757;
    // TODO: Check why MATLAB calculates b0, b1, b2 coefficients that slightly differ from 1.0 or 2.0 by about 0.01
    std::vector<double> expected_coeffs = {
            1.0, 1.0, 0.0, 1.0, 0.414213562373084, 0.0,
            1.0, 2.0, 1.0, 1.0, 0.863862810051667, 0.221686502004782,
            1.0, 2.0, 1.0, 1.0, 0.981497128142559, 0.388046550049447,
            1.0, 2.0, 1.0, 1.0, 1.22194528464162, 0.728091594018049,
    };
    normalize_coeffs(expected_coeffs);
    test_coeffs(butter_lpf, GAIN_EXPECTED, expected_coeffs, TOL);
}

TEST(ButterTest, ButterLPFQ31ImpulseResponse) {
    IIRButter<7, q31_t> butter_lpf(0.75);

    const std::vector<double> expected = {
            0.161071111239757, 0.566725679142607, 0.489525669314195, -0.192071409434096,
            -0.155083424057316, 0.219030913308335, -0.0872329542708532, -0.0648084553449979,
            0.128609514340678, -0.0937359884397627, 0.0120117665505027, 0.0552459536354045,
            -0.0741946034189801, 0.0477336794172441, -0.00248775660321518, -0.0324101273151006,
            0.0413620861236263, -0.0266064909895995, 0.00207667261047778, 0.0170202884196747,
    };

    test_impulse_response(butter_lpf, expected, TOL);
}
