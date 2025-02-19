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

TEST(EllipticTest, EllipticFunctions) {
    constexpr double TOL_EF = 1e-6;
    constexpr double ELLIPTIC_MODULUS = 0.5;

    const double K = calculate_elliptic_integral(ELLIPTIC_MODULUS);
    const double K_prime = calculate_elliptic_integral(get_complimentary(ELLIPTIC_MODULUS));
    const double R = K_prime / K;

    constexpr double ELLIPTIC_FUNCTION_ARG = 0.3;
    const Complex sn_val = sn(ELLIPTIC_FUNCTION_ARG, ELLIPTIC_MODULUS);
    const Complex cd_val = cd(ELLIPTIC_FUNCTION_ARG, ELLIPTIC_MODULUS);

    ASSERT_NEAR(sn_val.real(), 0.479967620395883, TOL_EF);
    ASSERT_NEAR(cd_val.real(), 0.903694981323325, TOL_EF);

    constexpr double ARC_ARG = 0.9;
    const Complex asn_val = asn(ARC_ARG, ELLIPTIC_MODULUS, R);
    const Complex acd_val = 1.0 - asn_val;

    ASSERT_NEAR(asn_val.real(), 0.694315905982350, TOL_EF);
    ASSERT_NEAR(acd_val.real(), 0.305684094017650, TOL_EF);

    const Complex sn_recalc = sn(asn_val, ELLIPTIC_MODULUS);
    const Complex cd_recalc = cd(acd_val, ELLIPTIC_MODULUS);

    ASSERT_NEAR(ARC_ARG, sn_recalc.real(), TOL_EF);
    ASSERT_NEAR(ARC_ARG, cd_recalc.real(), TOL_EF);
}
