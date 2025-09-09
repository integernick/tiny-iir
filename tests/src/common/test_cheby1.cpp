#include "test_utils.h"

#include <cheby1/IIRCheby1.h>

using namespace tiny_iir;

constexpr double TOL = TOL_DOUBLE;

TEST(Cheby1Test, Cheby1LPFDoubleCoeffs) {
    IIRCheby1<1, double> cheby1_lpf(0.75, 0.1);
    double GAIN_EXPECTED = 0.940541375070859;
    std::vector<double> expected_coeffs = {
            1.0, 1.0, 0.0, 1.0, 0.881082750141717, 0.0
    };
    test_coeffs(cheby1_lpf, GAIN_EXPECTED, expected_coeffs, TOL);
    
    IIRCheby1<4, double> cheby1_lpf_order_4(0.5, 0.3);
    GAIN_EXPECTED = 0.0765461714674103;
    expected_coeffs = {
            1, 2.0, 1.0, 1.0, -0.47674363011616, 0.186488199316368,
            1, 2.0, 1.0, 1.0, 0.105420579500914, 0.680825869653555
    };

    test_coeffs(cheby1_lpf_order_4, GAIN_EXPECTED, expected_coeffs, TOL);
}

TEST(Cheby1Test, Cheby1LPFDoubleImpulseResponse) {
    IIRCheby1<4, double> cheby1_lpf(0.5, 0.3);

    std::vector<double> expected_impulse_response = {
            0.0765461714674103, 0.334608043771838, 0.520982239568247, 0.249584500664994,
            -0.164138919075087, -0.148499040135320, 0.0889256922160354, 0.0726141799978015,
            -0.0701341279408741, -0.0394027086405953, 0.0535232951121069, 0.0214638234950325,
            -0.0388715022959289, -0.0106479082515899, 0.0275554646590436, 0.00435404845076304,
            -0.0192089830548209, -0.00093611102586206, 0.0131762335480102, -0.000752519593701615
    };

    test_impulse_response(cheby1_lpf, expected_impulse_response, TOL);

    cheby1_lpf.configure(0.7, 0.1);

    cheby1_lpf.reset_state(); // Resets filter state
    expected_impulse_response = {
            0.28752684851852, 0.682626808501851, 0.210378973651761, -0.300402215664622,
            0.0993629445718004, 0.0802308147374734, -0.139353842738623, 0.0951387939121066,
            -0.0102022082328173, -0.0547428256344923, 0.0697204148176845, -0.0413029677973009,
            -0.00192177659181128, 0.0315200844145639, -0.0349262988458463, 0.0178340885715921,
            0.00426928040051416, -0.0175849654811976, 0.0172198156255823, -0.00734039967472663,
    };

    test_impulse_response(cheby1_lpf, expected_impulse_response, TOL);
}

TEST(Cheby1Test, Cheby1BPFDoubleCoeffs) {
    IIRCheby1<4, double, FilterPassType::BAND_PASS> cheby1_bpf(0.7, 0.8, 0.1);

    constexpr double GAIN_EXPECTED = 0.000378930474568963;
    /** @note MATLAB sorts zeros so that positive ones come first */
    const std::vector<double> expected_coeffs = {
            1.0, 2.0, 1.0, 1, 1.19127613367568, 0.804935766025284,
            1.0, -2.0, 1.0, 1, 1.39810655164099, 0.830114885398173,
            1.0, 2.0, 1.0, 1, 1.08715913015724, 0.908892249014297,
            1.0, -2.0, 1.0, 1, 1.58373195838602, 0.93540445652867,
    };

    test_coeffs(cheby1_bpf, GAIN_EXPECTED, expected_coeffs, TOL);
}

TEST(Cheby1Test, Cheby1BPFImpulseResponse) {
    IIRCheby1<7, double, FilterPassType::BAND_PASS> cheby1_bpf(0.4, 0.6, 0.1);

    const std::vector<double> expected = {
            1.92703455044235e-05, 5.77648294330306e-20, -0.000237524848074349, -8.08921536809345e-19,
            0.00142415810127082, 5.52827323756361e-18, -0.00557286413357188, -2.46093089696208e-17,
            0.0161066755162555, 8.01547878464967e-17, -0.0368520506848039, -2.02498489925638e-16,
            0.0697297871057752, 4.07077099603553e-16, -0.112254446907077, -6.45945347751709e-16,
            0.15646607297403, 7.45022315216523e-16, -0.190393482833423, -3.58656327037608e-16
    };

    test_impulse_response(cheby1_bpf, expected, TOL);
}

TEST(Cheby1Test, Cheby1BSFImpulseResponse) {
    IIRCheby1<5, double, FilterPassType::BAND_STOP> cheby1_bsf(0.35, 0.45, 0.05);

    const std::vector<double> expected = {
            0.5962523346637, -0.191490997136767, 0.524010068943812, 0.332580223670483,
            -0.108109115141162, -0.171196354642223, -0.0197928045064232, 0.00495260703021452,
            -0.0368951567164588, 0.0319512096022865, 0.112772045022167, 0.0327751255331205,
            -0.112489728314173, -0.103498819058062, 0.0445415424213215, 0.114547592844524,
            0.0269861868090479, -0.0694601980209971, -0.0519178637498459, 0.0160700533743126
    };

    test_impulse_response(cheby1_bsf, expected, TOL);
}
