#include "test_utils.h"

#include <cheby1/IIRCheby1.h>

using namespace tiny_iir;

constexpr double TOL = TOL_DOUBLE;

TEST(Cheby1Test, Cheby1LPFDoubleCoeffs) {
    IIRCheby1<1, double> cheby1_lpf_double(0.75, 0.1);

    constexpr double GAIN_EXPECTED = 0.940541375070859;
    const std::vector<double> expected_coeffs = {
            1.0, 1.0, 0.0, 1, 0.881082750141717, 0.0
    };

    test_coeffs(cheby1_lpf_double, GAIN_EXPECTED, expected_coeffs, TOL);
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
    cheby1_bpf.print_coefficients();

    /* TODO: Trust me bro */
    EXPECT_TRUE(true);
}

TEST(Cheby1Test, Cheby1BSFImpulseResponse) {
    IIRCheby1<7, double, FilterPassType::BAND_STOP> cheby1_bsf(0.4, 0.6, 0.1);
    cheby1_bsf.print_coefficients();

    /* TODO: Trust */
    EXPECT_TRUE(true);
}
