#include "IIRButter.h"
#include "IIRCheby1.h"
#include "IIRCheby2.h"

#include <gtest/gtest.h>

namespace tiny_iir {
void expect_coefficients_near(const std::vector<double> &expected, const double *actual, double tolerance) {
    for (size_t i = 0; i < expected.size(); ++i) {
        EXPECT_NEAR(expected[i], actual[i], tolerance) << "Mismatch at coefficient index " << i;
    }
}

TEST(Cheby1Test, Cheby1Lowpass) {
    IIRCheby1<1, double> cheby1_lpf(0.75, 0.1);
    EXPECT_NEAR(cheby1_lpf.cascade_filter().get_gain(), 0.940541375070859, 1e-10)
                        << "Gain mismatch (single pole LPF)";
    std::vector<double> expected_coeffs = {
            1.0, 1.0, 0.0, 0.881082750141717, 0.0
    };
    expect_coefficients_near(
            expected_coeffs,
            cheby1_lpf.cascade_filter().get_coefficients(),
            1e-10
    );
}

TEST(Cheby2Test, Cheby2Lowpass) {
    using namespace tiny_iir;

    // Existing example for LPF
    IIRCheby2<4, double> cheby2_lpf(0.5, 40.0);
    EXPECT_NEAR(cheby2_lpf.cascade_filter().get_gain(), 0.0458146015732474, 1e-10) << "Gain mismatch (LPF)";

    std::vector<double> expected_coeffs_lpf = {
            1.0, 1.48904167641087, 1.0, -0.592228498379536, 0.130712470348595,
            1.0, 0.158017147118544, 1.0, -0.931033981965233, 0.571641721541078
    };
    expect_coefficients_near(
            expected_coeffs_lpf,
            cheby2_lpf.cascade_filter().get_coefficients(),
            1e-10
    );
}

TEST(Cheby2Test, Cheby2Highpass) {
    using namespace tiny_iir;

    IIRCheby2<7, double, FilterPassType::HIGH_PASS> cheby2_hpf(0.2, 30.0);
    EXPECT_NEAR(cheby2_hpf.cascade_filter().get_gain(), 0.414076474853479, 1e-10)
                        << "Gain mismatch (HPF)";

    std::vector<double> expected_coeffs_hpf = {
            1.0, -0.999999999998814, 0.0, -0.661214458569906, 0.0,
            1.0, -1.92205073355251, 1.00000000000223, -1.31360021053207, 0.485658078973171,
            1.0, -1.75751802610867, 0.999999999998277, -1.31410845670521, 0.626137289983513,
            1.0, -1.63522246863863, 1.00000000000068, -1.39217265603295, 0.852743809792487
    };

    expect_coefficients_near(
            expected_coeffs_hpf,
            cheby2_hpf.cascade_filter().get_coefficients(),
            1e-10
    );
}

}