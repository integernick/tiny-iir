#include "test_utils.h"
#include "test_utils_fixed_point.h"

#include "cheby2/IIRCheby2.h"

using namespace tiny_iir;

constexpr double TOL = TOL_Q31;

TEST(Cheby2Test, Cheby2LPFQ31Coeffs) {
    IIRCheby2<4, q31_t> cheby2_lpf(0.5, 40.0);

    constexpr double GAIN_EXPECTED = 0.0458146015732474;
    std::vector<double> expected_coeffs = {
            1.0, 1.48904167641087, 1.0, 1, -0.592228498379536, 0.130712470348595,
            1.0, 0.158017147118544, 1.0, 1, -0.931033981965233, 0.571641721541078,
    };
    normalize_coeffs(expected_coeffs);

    test_coeffs(cheby2_lpf, GAIN_EXPECTED, expected_coeffs, TOL);
}

TEST(Cheby2Test, Cheby2HPFQ31Coeffs) {
    IIRCheby2<7, q31_t, FilterPassType::HIGH_PASS> cheby2_hpf(0.2, 30.0);

    constexpr double GAIN_EXPECTED = 0.414076474853479;
    std::vector<double> expected_coeffs = {
            1.0, -0.999999999998814, 0.0, 1, -0.661214458569906, 0.0,
            1.0, -1.92205073355251, 1.00000000000223, 1, -1.31360021053207, 0.485658078973171,
            1.0, -1.75751802610867, 0.999999999998277, 1, -1.31410845670521, 0.626137289983513,
            1.0, -1.63522246863863, 1.00000000000068, 1, -1.39217265603295, 0.852743809792487,
    };
    normalize_coeffs(expected_coeffs);

    test_coeffs(cheby2_hpf, GAIN_EXPECTED, expected_coeffs, TOL_DOUBLE);
}
