#include "test_utils.h"
#include "test_utils_fixed_point.h"

#include <cheby1/IIRCheby1.h>

using namespace tiny_iir;

constexpr double TOL = TOL_Q31;

TEST(Cheby1Test, Cheby1LPFQ31Coeffs) {
    IIRCheby1<1, q31_t> cheby1_lpf_q31(0.75, 0.1);

    constexpr double GAIN_EXPECTED = 0.940541375070859;
    std::vector<double> expected_coeffs = {
            1.0, 1.0, 0.0, -0.881082750141717, 0.0
    };
    normalize_coeffs(expected_coeffs);

    test_coeffs(cheby1_lpf_q31, GAIN_EXPECTED, expected_coeffs, TOL);
}

TEST(Cheby1Test, Cheby1LPFQ31ImpulseResponse) {
    IIRCheby1<4, q31_t> cheby1_lpf(0.5, 0.3);

    std::vector<double> expected = {
            0.0765461714674103, 0.334608043771838, 0.520982239568247, 0.249584500664994,
            -0.164138919075087, -0.148499040135320, 0.0889256922160354, 0.0726141799978015,
            -0.0701341279408741, -0.0394027086405953, 0.0535232951121069, 0.0214638234950325,
            -0.0388715022959289, -0.0106479082515899, 0.0275554646590436, 0.00435404845076304,
            -0.0192089830548209, -0.00093611102586206, 0.0131762335480102, -0.000752519593701615
    };

    test_impulse_response(cheby1_lpf, expected, TOL);
    cheby1_lpf.reset();
    test_impulse_response_batch(cheby1_lpf, expected, TOL);
}