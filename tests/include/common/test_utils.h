#pragma once

#include <type_utils.h>
#include <common/CascadeFilter.h>

#include <gtest/gtest.h>
#include <type_traits>
#include <cmath>

namespace tiny_iir {
static constexpr double TOL_DOUBLE = 2e-9;

namespace {

template<class FILTER, typename DESIGN_T>
void test_biquad_coefficients(FILTER &filter, uint32_t biquad_idx, const BiquadCoefficients<DESIGN_T> &coeffs_expected,
                              double tolerance) {
    BiquadCoefficients<DESIGN_T> coeffs_actual = filter.get_biquad_coefficients(biquad_idx);

    EXPECT_NEAR(coeffs_expected.b0, coeffs_actual.b0, tolerance) << "Biquad coefficients mismatch";
    EXPECT_NEAR(coeffs_expected.b1, coeffs_actual.b1, tolerance) << "Biquad coefficients mismatch";
    EXPECT_NEAR(coeffs_expected.b2, coeffs_actual.b2, tolerance) << "Biquad coefficients mismatch";
    EXPECT_NEAR(coeffs_expected.a1, coeffs_actual.a1, tolerance) << "Biquad coefficients mismatch";
    EXPECT_NEAR(coeffs_expected.a2, coeffs_actual.a2, tolerance) << "Biquad coefficients mismatch";
}

} // namespace

template<class FILTER>
void test_coeffs(FILTER &filter, const double expected_gain, const std::vector<double> &expected_coeffs,
                 const double tolerance) {
    using ValueType = typename FILTER::ValueType;

    // arm_biquad_cascade_df2T_instance_f32/f64, arm_biquad_casd_df1_inst_q31 expect [b0, b1, b2, a1, a2]
    // arm_biquad_casd_df1_inst_q15 expects [b0, 0, b1, b2, a1, a2]
    constexpr uint32_t COEFFICIENTS_PER_BIQUAD_BLOCK_WITH_A0 = std::is_same_v<ValueType, q15_t>
                                                               ? coeffs_per_stage<ValueType>::value
                                                               : coeffs_per_stage<ValueType>::value + 1;

    const uint32_t num_of_biquad_blocks = (FILTER::REAL_ORDER + 1) / 2;
    const uint32_t num_of_coefficients = num_of_biquad_blocks * coeffs_per_stage<ValueType>::value;
    const uint32_t num_of_blocks_expected = expected_coeffs.size() / COEFFICIENTS_PER_BIQUAD_BLOCK_WITH_A0;
    const uint32_t num_of_coefficients_expected = std::is_same_v<ValueType, q15_t>
                                                  ? expected_coeffs.size()
                                                  : expected_coeffs.size() - num_of_blocks_expected;

    EXPECT_EQ(num_of_coefficients, num_of_coefficients_expected)
                        << "Expected coefficients size must be " << num_of_coefficients;
    EXPECT_EQ(num_of_biquad_blocks, num_of_blocks_expected)
                        << "Number of blocks mismatch";

    const ValueType gain = filter.get_gain();
    double gain_as_double;
    to_double(&gain, &gain_as_double, 1);

    EXPECT_NEAR(gain_as_double, expected_gain, tolerance) << "Filter gain mismatch";

    // Normalize coefficients if necessary
    for (uint32_t i = 0; i < num_of_biquad_blocks; ++i) {
        std::vector<double> expected_biquad_coeffs{expected_coeffs.begin() + i * COEFFICIENTS_PER_BIQUAD_BLOCK_WITH_A0,
                                                   expected_coeffs.begin() +
                                                   (i + 1) * COEFFICIENTS_PER_BIQUAD_BLOCK_WITH_A0};
        // a0 is always 1.0 so it is not stored in the actual coefficients
        constexpr uint32_t IDX_A0 = 3;
        expected_biquad_coeffs.erase(expected_biquad_coeffs.begin() + IDX_A0);

        BiquadCoefficients coeffs_expected{
                .b0 = expected_biquad_coeffs[0],
                .b1 = expected_biquad_coeffs[1],
                .b2 = expected_biquad_coeffs[2],
                .a1 = expected_biquad_coeffs[3],
                .a2 = expected_biquad_coeffs[4],
        };

        test_biquad_coefficients(filter, i, coeffs_expected, tolerance);
    }
}

template<class FILTER>
void test_response(FILTER &filter, const std::vector<double> &input_signal,
                   const std::vector<double> &expected_response, const double tolerance) {
    const uint32_t response_size = expected_response.size();
    std::vector<double> response(response_size);

    filter.reset_state();
    for (int i = 0; i < input_signal.size(); ++i) {
        response[i] = filter.process(input_signal[i]);
    }

    for (uint32_t i = 0; i < response_size; ++i) {
        EXPECT_NEAR(response[i], expected_response[i], tolerance) << "Mismatch at sample index " << i;
    }

    // Test batch processing
    filter.reset_state();
    const double last_output = filter.process(input_signal.data(), input_signal.size());
    EXPECT_NEAR(last_output, expected_response.back(), tolerance) << "Mismatch in final returned sample";
}

template<class FILTER>
void test_impulse_response(FILTER &filter, const std::vector<double> &expected_response, const double tolerance) {
    const uint32_t response_size = expected_response.size();
    std::vector<double> impulse(response_size);
    impulse[0] = 1.0;

    test_response(filter, impulse, expected_response, tolerance);
}

template<class FILTER>
void test_step_response(FILTER &filter, const std::vector<double> &expected_response, const double tolerance) {
    const uint32_t response_size = expected_response.size();
    std::vector<double> step(response_size);
    std::fill(step.begin(), step.end(), 1.0);

    test_response(filter, step, expected_response, tolerance);
}

}