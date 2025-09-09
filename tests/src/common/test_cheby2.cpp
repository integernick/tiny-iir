#include "test_utils.h"

#include "cheby2/IIRCheby2.h"

using namespace tiny_iir;

constexpr double TOL = TOL_DOUBLE;

TEST(Cheby2Test, Cheby2LPFDoubleCoeffs) {
    IIRCheby2<4, double> cheby2_lpf(0.5, 40.0);

    constexpr double GAIN_EXPECTED = 0.0458146015732474;
    const std::vector<double> expected_coeffs = {
            1.0, 1.48904167641087, 1.0, 1, -0.592228498379536, 0.130712470348595,
            1.0, 0.158017147118544, 1.0, 1, -0.931033981965233, 0.571641721541078,
    };

    test_coeffs(cheby2_lpf, GAIN_EXPECTED, expected_coeffs, TOL);
}

TEST(Cheby2Test, Cheby2HPFDoubleCoeffs) {
    IIRCheby2<7, double, FilterPassType::HIGH_PASS> cheby2_hpf(0.2, 30.0);
    EXPECT_NEAR(cheby2_hpf.get_gain(), 0.414076474853479, 1e-10) << "Gain mismatch (HPF)";

    constexpr double GAIN_EXPECTED = 0.414076474853479;
    const std::vector<double> expected_coeffs = {
            1.0, -0.999999999998814, 0.0, 1, -0.661214458569906, 0.0,
            1.0, -1.92205073355251, 1.00000000000223, 1, -1.31360021053207, 0.485658078973171,
            1.0, -1.75751802610867, 0.999999999998277, 1, -1.31410845670521, 0.626137289983513,
            1.0, -1.63522246863863, 1.00000000000068, 1, -1.39217265603295, 0.852743809792487,
    };

    test_coeffs(cheby2_hpf, GAIN_EXPECTED, expected_coeffs, TOL);
}

TEST(Cheby2Test, Cheby2BPFImpulseResponse) {
    IIRCheby2<7, double, FilterPassType::BAND_PASS> cheby2_bpf(0.05, 0.35, 30.0);

    const std::vector<double> expected = {
            0.0513292514145267, 0.0963269056660768, 0.122911666862606, 0.104800727679869,
            0.00027661820339131, -0.116167702046018, -0.220476511215389, -0.24640929389168,
            -0.167438648466226, -0.0330693742280879, 0.0785403657897334, 0.117976891120525,
            0.0876663814909757, 0.0304242353789907, -0.00280252898017574, 0.0123385068009629,
            0.0584958195244826, 0.0933015029883063, 0.0854650264815308, 0.0401324286533012
    };

    test_impulse_response(cheby2_bpf, expected, TOL);
}

TEST(Cheby2Test, Cheby2BSFImpulseResponse) {
    IIRCheby2<7, double, FilterPassType::BAND_STOP> cheby2_bsf(0.05, 0.35, 30.0);

    const std::vector<double> expected = {
            0.268494049405138, -0.539362745210064, 0.381062283143844, 0.130409607636981,
            -0.164043905376265, -0.0883830841165492, 0.107484437394004, 0.128637930747905,
            0.00275375027782676, -0.066287926668191, -0.00128535400870133, 0.0879712870934486,
            0.086093165397145, 0.0154262486583965, -0.0209432585406577, 0.0176602949773436,
            0.0715483225421359, 0.07264425625174, 0.0297999546789769, 0.00498295964691192
    };

    test_impulse_response(cheby2_bsf, expected, TOL);
}
