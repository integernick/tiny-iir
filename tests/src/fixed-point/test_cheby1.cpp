#include "test_utils.h"
#include "test_utils_fixed_point.h"

#include <cheby1/IIRCheby1.h>

using namespace tiny_iir;

constexpr double TOL = TOL_Q31;

TEST(Cheby1Test, Cheby1LPFQ31Coeffs) {
    IIRCheby1<1, q31_t> cheby1_lpf_q31(0.75, 0.1);

    constexpr double GAIN_EXPECTED = 0.940541375070859;
    std::vector<double> expected_coeffs = {
            1.0, 1.0, 0.0, 1, 0.881082750141717, 0.0
    };
    normalize_coeffs(expected_coeffs);

    test_coeffs(cheby1_lpf_q31, GAIN_EXPECTED, expected_coeffs, TOL);
}

TEST(Cheby1Test, Cheby1LPFQ31ImpulseResponse) {
    IIRCheby1<4, q31_t> cheby1_lpf(0.5, 0.3);

    const std::vector<double> expected = {
            0.0765461714674103, 0.334608043771838, 0.520982239568247, 0.249584500664994, -0.164138919075087,
            -0.148499040135320, 0.0889256922160354, 0.0726141799978015, -0.0701341279408741, -0.0394027086405953,
            0.0535232951121069, 0.0214638234950325, -0.0388715022959289, -0.0106479082515899, 0.0275554646590436,
            0.00435404845076304, -0.0192089830548209, -0.00093611102586206, 0.0131762335480102, -0.000752519593701615
    };

    test_impulse_response(cheby1_lpf, expected, TOL);
}

TEST(Cheby1Test, Cheby1BPFQ31Coeffs) {
    IIRCheby1<4, q31_t, FilterPassType::BAND_PASS> cheby1_bpf(0.7, 0.8, 0.1);

    constexpr double GAIN_EXPECTED = 0.000378930474568963;
    /** @note MATLAB sorts zeros so that positive ones come first */
    std::vector<double> expected_coeffs = {
            1.0, 2.0, 1.0, 1, 1.19127613367568, 0.804935766025284,
            1.0, -2.0, 1.0, 1, 1.39810655164099, 0.830114885398173,
            1.0, 2.0, 1.0, 1, 1.08715913015724, 0.908892249014297,
            1.0, -2.0, 1.0, 1, 1.58373195838602, 0.93540445652867,
    };
    normalize_coeffs(expected_coeffs);

    test_coeffs(cheby1_bpf, GAIN_EXPECTED, expected_coeffs, TOL);
}

TEST(Cheby1Test, Cheby1BPFQ31ImpulseResponse) {
    IIRCheby1<7, q31_t, FilterPassType::BAND_PASS> cheby1_bpf(0.4, 0.6, 0.1);

    const std::vector<double> expected = {
            1.92703455044235e-05, 5.77648294330306e-20, -0.000237524848074349, -8.08921536809345e-19,
            0.00142415810127082, 5.52827323756361e-18, -0.00557286413357188, -2.46093089696208e-17,
            0.0161066755162555, 8.01547878464967e-17, -0.0368520506848039, -2.02498489925638e-16,
            0.0697297871057752, 4.07077099603553e-16, -0.112254446907077, -6.45945347751709e-16,
            0.15646607297403, 7.45022315216523e-16, -0.190393482833423, -3.58656327037608e-16
    };

    test_impulse_response(cheby1_bpf, expected, 1e-4);
}

TEST(Cheby1Test, Cheby1BSFQ31ImpulseResponse) {
    IIRCheby1<5, q31_t, FilterPassType::BAND_STOP> cheby1_bsf(0.35, 0.45, 0.05);

    const std::vector<double> expected = {
            0.5962523346637, -0.191490997136767, 0.524010068943812, 0.332580223670483,
            -0.108109115141162, -0.171196354642223, -0.0197928045064232, 0.00495260703021452,
            -0.0368951567164588, 0.0319512096022865, 0.112772045022167, 0.0327751255331205,
            -0.112489728314173, -0.103498819058062, 0.0445415424213215, 0.114547592844524,
            0.0269861868090479, -0.0694601980209971, -0.0519178637498459, 0.0160700533743126
    };

    test_impulse_response(cheby1_bsf, expected, 1e-4);
}