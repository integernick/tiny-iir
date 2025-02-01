#include "test_utils.h"

#include "cheby2/IIRCheby2.h"

using namespace tiny_iir;

constexpr double TOL = TOL_DOUBLE;

TEST(Cheby2Test, Cheby2LPFDoubleCoeffs) {
    IIRCheby2<4, double> cheby2_lpf(0.5, 40.0);

    constexpr double GAIN_EXPECTED = 0.0458146015732474;
    const std::vector<double> expected_coeffs = {
            1.0, 1.48904167641087, 1.0, 0.592228498379536, -0.130712470348595,
            1.0, 0.158017147118544, 1.0, 0.931033981965233, -0.571641721541078
    };

    test_coeffs(cheby2_lpf, GAIN_EXPECTED, expected_coeffs, TOL);
}

TEST(Cheby2Test, Cheby2HPFDoubleCoeffs) {
    IIRCheby2<7, double, FilterPassType::HIGH_PASS> cheby2_hpf(0.2, 30.0);
    EXPECT_NEAR(cheby2_hpf.cascade_filter().get_gain(), 0.414076474853479, 1e-10)
                        << "Gain mismatch (HPF)";

    constexpr double GAIN_EXPECTED = 0.414076474853479;
    const std::vector<double> expected_coeffs = {
            1.0, -0.999999999998814, 0.0, 0.661214458569906, 0.0,
            1.0, -1.92205073355251, 1.00000000000223, 1.31360021053207, -0.485658078973171,
            1.0, -1.75751802610867, 0.999999999998277, 1.31410845670521, -0.626137289983513,
            1.0, -1.63522246863863, 1.00000000000068, 1.39217265603295, -0.852743809792487
    };

    test_coeffs(cheby2_hpf, GAIN_EXPECTED, expected_coeffs, TOL);
}
