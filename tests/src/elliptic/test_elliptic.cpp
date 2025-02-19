#include "test_utils.h"

#include <elliptic/IIRElliptic.h>

using namespace tiny_iir;

constexpr double TOL = TOL_DOUBLE;

TEST(EllipticTest, EllipticLPFCoeffs) {
    IIRElliptic<5, double, FilterPassType::HIGH_PASS> lpf_elliptic(0.4, 0.05, 40.0);
    lpf_elliptic.print_coefficients();
}

TEST(EllipticTest, EllipticSolver) {
    constexpr double TOL_EF = 1e-6;
    constexpr double ELLIPTIC_MODULUS = 0.5;

    EllipticSolver solver;
    solver.init(ELLIPTIC_MODULUS);

    constexpr double ELLIPTIC_FUNCTION_ARG = 0.3;
    const double sn = solver.sn(ELLIPTIC_FUNCTION_ARG);
    const double cd = solver.cd(ELLIPTIC_FUNCTION_ARG);

    ASSERT_NEAR(sn, 0.479967620395883, TOL_EF);
    ASSERT_NEAR(cd, 0.903694981323325, TOL_EF);

    constexpr double ARC_ARG = 0.9;
    const double asn = solver.asn(ARC_ARG);
    const double acd = solver.acd(ARC_ARG);

    ASSERT_NEAR(asn, 0.694315905982350, TOL_EF);
    ASSERT_NEAR(acd, 0.305684094017650, TOL_EF);

    const double sn_recalc = solver.sn(asn);
    const double cd_recalc = solver.cd(acd);

    ASSERT_NEAR(ARC_ARG, sn_recalc, TOL_EF);
    ASSERT_NEAR(ARC_ARG, cd_recalc, TOL_EF);
}
