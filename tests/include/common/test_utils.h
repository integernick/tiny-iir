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

    constexpr size_t NUM_OF_BIQUAD_BLOCKS = CascadeFilter<FILTER::ORDER, ValueType>::NUMBER_OF_BIQUAD_BLOCKS;
    constexpr size_t NUM_OF_COEFFICIENTS = CascadeFilter<FILTER::ORDER, ValueType>::NUMBER_OF_COEFFICIENTS;
    constexpr size_t COEFFICIENTS_PER_BIQUAD_BLOCK = CascadeFilter<FILTER::ORDER, ValueType>::COEFFICIENTS_PER_BIQUAD_BLOCK;

    EXPECT_EQ(expected_coeffs.size(), NUM_OF_COEFFICIENTS)
                        << "Expected coefficients size must be " << NUM_OF_COEFFICIENTS;
    EXPECT_EQ(filter.cascade_filter().get_number_of_blocks(), expected_coeffs.size() / COEFFICIENTS_PER_BIQUAD_BLOCK)
                        << "Number of blocks mismatch";

    const ValueType gain = filter.cascade_filter().get_gain();
    double gain_as_double;
    to_double(&gain, &gain_as_double, 1);

    EXPECT_NEAR(gain_as_double, expected_gain, tolerance) << "Filter gain mismatch";

    // Normalize coefficients if necessary
    for (size_t i = 0; i < NUM_OF_BIQUAD_BLOCKS; ++i) {
        std::vector<double> expected_biquad_coeffs{expected_coeffs.begin() + i * COEFFICIENTS_PER_BIQUAD_BLOCK,
                                                   expected_coeffs.begin() + (i + 1) * COEFFICIENTS_PER_BIQUAD_BLOCK};

        for (size_t j = 0; j < COEFFICIENTS_PER_BIQUAD_BLOCK; ++j) {
            const size_t coeff_idx = i * COEFFICIENTS_PER_BIQUAD_BLOCK + j;
            double coeff_as_double;
            to_double(&filter.cascade_filter().get_coefficients()[coeff_idx], &coeff_as_double, 1);
            EXPECT_NEAR(expected_biquad_coeffs[j], coeff_as_double, tolerance) << "Mismatch at coefficient index " << i;
        }
    }
}

template<class FILTER>
void test_impulse_response(FILTER &filter, const std::vector<double> &expected_response, const double tolerance) {
    using ValueType = typename FILTER::ValueType;

    const size_t response_size = expected_response.size();
    std::vector<double> impulse(response_size);
    impulse[0] = 1.0;

    std::vector<ValueType> impulse_response(response_size);

    for (int i = 0; i < impulse.size(); ++i) {
        impulse_response[i] = filter.process(impulse[i]);
    }

    for (size_t i = 0; i < response_size; ++i) {
        double sample_as_double;
        to_double(&impulse_response[i], &sample_as_double, 1);
        EXPECT_NEAR(sample_as_double, expected_response[i], tolerance) << "Mismatch at sample index " << i;
    }
}

template<class FILTER>
void test_impulse_response_batch(FILTER &filter, const std::vector<double> &expected_response, const double tolerance) {
    const size_t response_size = expected_response.size();
    std::vector<double> impulse(response_size);
    impulse[0] = 1.0;

    const auto last_output = filter.process(impulse.data(), response_size);
    double last_output_as_double;
    to_double(&last_output, &last_output_as_double, 1);
    EXPECT_NEAR(last_output_as_double, expected_response.back(), tolerance) << "Mismatch in final returned sample";
}

}