#include "test_utils.h"

#include <butter/IIRButter.h>

using namespace tiny_iir;

constexpr double TOL = TOL_DOUBLE;

TEST(ButterTest, ButterLPFDoubleCoeffs) {
    IIRButter<7, double> butter_lpf(0.75);

    constexpr double GAIN_EXPECTED = 0.161071111239757;
    // TODO: Check why MATLAB calculates b0, b1, b2 coefficients that slightly differ from 1.0 or 2.0 by about 0.01
    const std::vector<double> expected_coeffs = {
            1.0, 1.0, 0.0, 1.0, 0.414213562373084, 0.0,
            1.0, 2.0, 1.0, 1.0, 0.863862810051667, 0.221686502004782,
            1.0, 2.0, 1.0, 1.0, 0.981497128142559, 0.388046550049447,
            1.0, 2.0, 1.0, 1.0, 1.22194528464162, 0.728091594018049,
    };

    test_coeffs(butter_lpf, GAIN_EXPECTED, expected_coeffs, TOL);
}

TEST(ButterTest, ButterLPFDoubleImpulseResponse) {
    IIRButter<7, double> butter_lpf(0.75);

    const std::vector<double> expected = {
            0.161071111239757, 0.566725679142607, 0.489525669314195, -0.192071409434096,
            -0.155083424057316, 0.219030913308335, -0.0872329542708532, -0.0648084553449979,
            0.128609514340678, -0.0937359884397627, 0.0120117665505027, 0.0552459536354045,
            -0.0741946034189801, 0.0477336794172441, -0.00248775660321518, -0.0324101273151006,
            0.0413620861236263, -0.0266064909895995, 0.00207667261047778, 0.0170202884196747,
    };

    test_impulse_response(butter_lpf, expected, TOL);
}

TEST(ButterTest, ButterBPFImpulseResponse) {
    IIRButter<4, double, FilterPassType::BandPass> butter_bpf(0.7, 0.8);

    const std::vector<double> expected = {
            0.0004165992044066, -0.00214164053195773, 0.00387518458087297, -0.000533796684804137,
            -0.0102354323971724, 0.0200022228497397, -0.0140686473219131, -0.0123543268463026,
            0.0432873295131373, -0.0506195705416532, 0.0180264930410516, 0.0402795170219338,
            -0.0845763656196623, 0.0777583143245449, -0.0155612122321824, -0.0657902132128914,
            0.112505913505717, -0.0909008480660411, 0.0120788251721443, 0.0747412773328065
    };

    test_impulse_response(butter_bpf, expected, TOL);
}

TEST(ButterTest, ButterBSFImpulseResponse) {
    IIRButter<4, double, FilterPassType::BandStop> butter_bsf(0.7, 0.8);

    const std::vector<double> expected = {
            0.662015837202617, 0.388337188268718, 0.105855702753534, -0.354102300634946,
            0.315115487960205, -0.131696469542741, -0.0218919652904209, 0.0640218474555734,
            -0.0273551091359787, -0.0069268835670151, -0.00938722354039277, 0.0588665355969819,
            -0.0875264038456233, 0.0592027602542121, 0.0130851922630142, -0.0801948810152771,
            0.0976008575545104, -0.0583155792486918, -0.00613558767237893, 0.0529505177717694
    };

    test_impulse_response(butter_bsf, expected, TOL);
}
