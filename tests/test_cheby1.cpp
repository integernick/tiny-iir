#include "IIRCheby1.h"
#include "test_utils.h"

using namespace tiny_iir;

static constexpr double TOL_Q31 = 1e-8;
static constexpr double TOL_DOUBLE = 1e-15;

TEST(Cheby1Test, Cheby1LowpassCoeffs) {
    IIRCheby1<1, double> cheby1_lpf_double(0.75, 0.1);
    IIRCheby1<1, q31_t> cheby1_lpf_q31(0.75, 0.1);

    const q31_t gain = cheby1_lpf_q31.cascade_filter().get_gain();
    double gain_q31_as_double;
    arm_q31_to_f64(&gain, &gain_q31_as_double, 1);

    constexpr double GAIN_EXPECTED = 0.940541375070859;
    EXPECT_NEAR(cheby1_lpf_double.cascade_filter().get_gain(), GAIN_EXPECTED, TOL_DOUBLE)
                        << "Gain (double) mismatch (single pole LPF)";
    EXPECT_NEAR(gain_q31_as_double, GAIN_EXPECTED, TOL_Q31)
                        << "Gain (q31) mismatch (single pole LPF)";

    const std::vector<double> expected_coeffs = {
            1.0, 1.0, 0.0, -0.881082750141717, 0.0
    };

    for (size_t i = 0; i < expected_coeffs.size(); ++i) {
        EXPECT_NEAR(expected_coeffs[i], cheby1_lpf_double.cascade_filter().get_coefficients()[i], TOL_DOUBLE)
                            << "Mismatch at coefficient index " << i;

        double actual_double;
        arm_q31_to_f64(&cheby1_lpf_q31.cascade_filter().get_coefficients()[i], &actual_double, 1);
        EXPECT_NEAR(expected_coeffs.at(i), actual_double, TOL_Q31) << "Mismatch at coefficient index " << i;
    }

}

TEST(Cheby1Test, Cheby1LowpassDoubleImpulseResponse) {
    IIRCheby1<4, double> cheby1_lpf_double(0.5, 0.3);
    IIRCheby1<4, q31_t> cheby1_lpf_q31(0.5, 0.3);

    constexpr size_t N = 20;
    std::array<double, 20> impulse{};
    impulse[0] = 1.0;

    std::vector<double> impulse_response_double(N);
    std::vector<q31_t> impulse_response_q31(N);

    for (int i = 0; i < impulse.size(); ++i) {
        impulse_response_double[i] = cheby1_lpf_double.process(impulse[i]);
        impulse_response_q31[i] = cheby1_lpf_q31.process(impulse[i]);
    }

    const std::vector<double> expected = {
            0.0765461714674103, 0.334608043771838, 0.520982239568247, 0.249584500664994,
            -0.164138919075087, -0.148499040135320, 0.0889256922160354, 0.0726141799978015,
            -0.0701341279408741, -0.0394027086405953, 0.0535232951121069, 0.0214638234950325,
            -0.0388715022959289, -0.0106479082515899, 0.0275554646590436, 0.00435404845076304,
            -0.0192089830548209, -0.00093611102586206, 0.0131762335480102, -0.000752519593701615
    };

    for (size_t i = 0; i < expected.size(); ++i) {
        EXPECT_NEAR(impulse_response_double[i], expected[i], TOL_DOUBLE) << "Mismatch at sample index " << i;
        double impulse_response_q31_as_double;
        arm_q31_to_f64(&impulse_response_q31[i], &impulse_response_q31_as_double, 1);
        EXPECT_NEAR(impulse_response_q31_as_double, expected[i], TOL_Q31) << "Mismatch at sample index " << i;
    }

    cheby1_lpf_double.reset();
    const double last_output_double = cheby1_lpf_double.process(impulse.data(), N);
    EXPECT_NEAR(last_output_double, expected.back(), TOL_DOUBLE) << "Mismatch in final returned sample";

    cheby1_lpf_q31.reset();
    const q31_t last_output_q31 = cheby1_lpf_q31.process(impulse.data(), N);
    double last_output_q31_as_double;
    arm_q31_to_f64(&last_output_q31, &last_output_q31_as_double, 1);
    EXPECT_NEAR(last_output_q31_as_double, expected.back(), TOL_Q31) << "Mismatch in final returned sample";
}
