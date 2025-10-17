#include "test_utils.h"
#include "test_utils_fixed_point.h"

#include <elliptic/IIRElliptic.h>

using namespace tiny_iir;

TEST(EllipticTest, EllipticLPFQ31Coeffs) {
    IIRElliptic<5, q31_t> lpf_elliptic(0.75, 0.05, 40.0);

    constexpr double GAIN_EXPECTED = 0.358219308431662;
    std::vector<double> expected_coeffs = {
            1.0, 1.0,              0,   1.0, 0.294668540839501, 0,
            1.0, 1.87575522417108, 1.0, 1.0, 0.992804485390727, 0.462892556367584,
            1.0, 1.73245210764594, 1.0, 1.0, 1.399394327975420, 0.860444958845612,
    };
    normalize_coeffs(expected_coeffs);

    test_coeffs(lpf_elliptic, GAIN_EXPECTED, expected_coeffs, TOL_Q31_GRADE2);
}

TEST(EllipticTest, EllipticLPFIQ31mpulseResponse) {
    IIRElliptic<5, q31_t> lpf_elliptic(0.75, 0.05, 40.0);

    const std::vector<double> expected = {
            0.358219308431662, 0.688241078042038, 0.0995525655880929, -0.27099728754872,
            0.153495101278209, 0.0220480430146786, -0.124648540267716, 0.122511675952959,
            -0.0511855851963433, -0.0308569973914753, 0.0781238804169214, -0.0750460151030332,
            0.0343210700065087, 0.0164223988543708, -0.05078383546014, 0.0552763038043805,
            -0.0328092958190739, -0.0017220189873963, 0.0303205522944549, -0.0405974147825837
    };

    test_impulse_response(lpf_elliptic, expected, TOL_Q31_GRADE2);
}

TEST(EllipticTest, EllipticLPFQ15Coeffs) {
    IIRElliptic<5, q15_t> lpf_elliptic(0.75, 0.05, 40.0);

    constexpr double GAIN_EXPECTED = 0.358219308431662;
    std::vector<double> expected_coeffs = {
            1.0, 1.0,              0,   1.0, 0.294668540839501, 0,
            1.0, 1.87575522417108, 1.0, 1.0, 0.992804485390727, 0.462892556367584,
            1.0, 1.73245210764594, 1.0, 1.0, 1.399394327975420, 0.860444958845612,
    };
    normalize_coeffs(expected_coeffs);

    test_coeffs(lpf_elliptic, GAIN_EXPECTED, expected_coeffs, TOL_Q15);
}

TEST(EllipticTest, EllipticLPFIQ15mpulseResponse) {
    IIRElliptic<5, q15_t> lpf_elliptic(0.75, 0.05, 40.0);

    const std::vector<double> expected = {
            0.358219308431662, 0.688241078042038, 0.0995525655880929, -0.27099728754872,
            0.153495101278209, 0.0220480430146786, -0.124648540267716, 0.122511675952959,
            -0.0511855851963433, -0.0308569973914753, 0.0781238804169214, -0.0750460151030332,
            0.0343210700065087, 0.0164223988543708, -0.05078383546014, 0.0552763038043805,
            -0.0328092958190739, -0.0017220189873963, 0.0303205522944549, -0.0405974147825837
    };

    test_impulse_response(lpf_elliptic, expected, TOL_Q15);
}
