#include "test_utils.h"

#include <elliptic/IIRElliptic.h>

using namespace tiny_iir;

constexpr double TOL = 1e-3; // Coefficients are not an exact match to MATLAB, but the frequency response is

TEST(EllipticTest, EllipticLPFCoeffs) {
    IIRElliptic<5, double> lpf_elliptic(0.75, 0.05, 40.0);

    constexpr double GAIN_EXPECTED = 0.358219308431662;
    std::vector<double> expected_coeffs = {
            1.0, 1.0,              0,   1.0, 0.294668540839501, 0,
            1.0, 1.87575522417108, 1.0, 1.0, 0.992804485390727, 0.462892556367584,
            1.0, 1.73245210764594, 1.0, 1.0, 1.399394327975420, 0.860444958845612,
    };

    test_coeffs(lpf_elliptic, GAIN_EXPECTED, expected_coeffs, TOL);
}

TEST(EllipticTest, EllipticLPFImpulseResponse) {
    IIRElliptic<5, double> lpf_elliptic(0.75, 0.05, 40.0);

    const std::vector<double> expected = {
            0.358219308431662, 0.688241078042038, 0.0995525655880929, -0.27099728754872,
            0.153495101278209, 0.0220480430146786, -0.124648540267716, 0.122511675952959,
            -0.0511855851963433, -0.0308569973914753, 0.0781238804169214, -0.0750460151030332,
            0.0343210700065087, 0.0164223988543708, -0.05078383546014, 0.0552763038043805,
            -0.0328092958190739, -0.0017220189873963, 0.0303205522944549, -0.0405974147825837
    };

    test_impulse_response(lpf_elliptic, expected, TOL);
}

TEST(EllipticTest, EllipticFunctions) {
    constexpr double TOL_EF = 1e-6;
    constexpr double ELLIPTIC_MODULUS = 0.5;

    const double K = calculate_elliptic_integral(ELLIPTIC_MODULUS);
    const double K_prime = calculate_elliptic_integral(get_complimentary(ELLIPTIC_MODULUS));
    const double R = K_prime / K;

    constexpr double ELLIPTIC_FUNCTION_ARG = 0.3;
    const Complex sn_val = sn<double>(ELLIPTIC_FUNCTION_ARG, ELLIPTIC_MODULUS);
    const Complex cd_val = cd<double>(ELLIPTIC_FUNCTION_ARG, ELLIPTIC_MODULUS);

    ASSERT_NEAR(sn_val.real(), 0.479967620395883, TOL_EF);
    ASSERT_NEAR(cd_val.real(), 0.903694981323325, TOL_EF);

    constexpr double ARC_ARG = 0.9;
    const Complex asn_val = asn<double>(ARC_ARG, ELLIPTIC_MODULUS, R);
    const Complex acd_val = 1.0 - asn_val;

    ASSERT_NEAR(asn_val.real(), 0.694315905982350, TOL_EF);
    ASSERT_NEAR(acd_val.real(), 0.305684094017650, TOL_EF);

    const Complex sn_recalc = sn(asn_val, ELLIPTIC_MODULUS);
    const Complex cd_recalc = cd(acd_val, ELLIPTIC_MODULUS);

    ASSERT_NEAR(ARC_ARG, sn_recalc.real(), TOL_EF);
    ASSERT_NEAR(ARC_ARG, cd_recalc.real(), TOL_EF);
}
