#pragma once

#include <utils.h>
#include <CascadeFilter.h>

#include <gtest/gtest.h>

namespace tiny_iir {

static constexpr double TOL_DOUBLE = 1e-9;

template<class FILTER>
void test_coeffs(FILTER &filter, const double expected_gain, const std::vector<double> &expected_coeffs,
                 const double tolerance) {
    using ValueType = typename FILTER::ValueType;

    constexpr size_t COEFFICIENTS_PER_BIQUAD_BLOCK = CascadeFilter<FILTER::ORDER, ValueType>::COEFFICIENTS_PER_BIQUAD_BLOCK;
    constexpr size_t COEFFICIENTS_PER_BIQUAD_BLOCK_WITH_A0 = COEFFICIENTS_PER_BIQUAD_BLOCK + 1;

    const size_t num_of_biquad_blocks = filter.get_number_of_blocks();
    const size_t num_of_coefficients = filter.get_number_of_blocks() * COEFFICIENTS_PER_BIQUAD_BLOCK;
    const size_t num_of_coefficients_expected = expected_coeffs.size() - num_of_biquad_blocks;
    const size_t num_of_blocks_expected = expected_coeffs.size() / COEFFICIENTS_PER_BIQUAD_BLOCK_WITH_A0;

    EXPECT_EQ(num_of_coefficients, num_of_coefficients_expected)
                        << "Expected coefficients size must be " << num_of_coefficients;
    EXPECT_EQ(num_of_biquad_blocks, num_of_blocks_expected)
                        << "Number of blocks mismatch";

    const ValueType gain = filter.get_gain();
    double gain_as_double;
    to_double(&gain, &gain_as_double, 1);

    EXPECT_NEAR(gain_as_double, expected_gain, tolerance) << "Filter gain mismatch";

    // Normalize coefficients if necessary
    for (size_t i = 0; i < num_of_biquad_blocks; ++i) {
        std::vector<double> expected_biquad_coeffs{expected_coeffs.begin() + i * COEFFICIENTS_PER_BIQUAD_BLOCK_WITH_A0,
                                                   expected_coeffs.begin() + (i + 1) * COEFFICIENTS_PER_BIQUAD_BLOCK_WITH_A0};
        // a0 is always 1.0 so it is not stored in the actual coefficients
        constexpr size_t IDX_A0 = 3;
        expected_biquad_coeffs.erase(expected_biquad_coeffs.begin() + IDX_A0);

        for (size_t j = 0; j < COEFFICIENTS_PER_BIQUAD_BLOCK; ++j) {
            const size_t coeff_idx = i * COEFFICIENTS_PER_BIQUAD_BLOCK + j;
            double coeff_expected = expected_biquad_coeffs[j];
            if (j > 2) {
                coeff_expected = -coeff_expected; // a0 and a1 are stored with an inverse sign compared to MATLAB
            }
            double coeff_as_double;
            to_double(&filter.get_coefficients()[coeff_idx], &coeff_as_double, 1);
            EXPECT_NEAR(coeff_expected, coeff_as_double, tolerance) << "Mismatch at coefficient index " << i;
        }
    }
}

template<class FILTER>
void test_response(FILTER &filter, const std::vector<double> &input_signal,
                   const std::vector<double> &expected_response, const double tolerance) {
    using ValueType = typename FILTER::ValueType;

    const size_t response_size = expected_response.size();
    std::vector<ValueType> response(response_size);

    for (int i = 0; i < input_signal.size(); ++i) {
        response[i] = filter.process(input_signal[i]);
    }

    for (size_t i = 0; i < response_size; ++i) {
        double sample_as_double;
        to_double(&response[i], &sample_as_double, 1);
        EXPECT_NEAR(sample_as_double, expected_response[i], tolerance) << "Mismatch at sample index " << i;
    }

    // Test batch processing
    filter.reset();
    const auto last_output = filter.process(input_signal.data(), input_signal.size());
    double last_output_as_double;
    to_double(&last_output, &last_output_as_double, 1);
    EXPECT_NEAR(last_output_as_double, expected_response.back(), tolerance) << "Mismatch in final returned sample";
}

template<class FILTER>
void test_impulse_response(FILTER &filter, const std::vector<double> &expected_response, const double tolerance) {
    using ValueType = typename FILTER::ValueType;

    const size_t response_size = expected_response.size();
    std::vector<double> impulse(response_size);
    impulse[0] = 1.0;

    test_response(filter, impulse, expected_response, tolerance);
}

template<class FILTER>
void test_step_response(FILTER &filter, const std::vector<double> &expected_response, const double tolerance) {
    using ValueType = typename FILTER::ValueType;

    const size_t response_size = expected_response.size();
    std::vector<double> step(response_size);
    std::fill(step.begin(), step.end(), 1.0);

    test_response(filter, step, expected_response, tolerance);
}

}