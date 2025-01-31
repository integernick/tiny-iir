#pragma once

#include <dsp/filtering_functions.h>
#include <string>

namespace tiny_iir {

template<typename T>
static void to_double(const T *x, double *x_double, size_t num_samples) {
    if constexpr (std::is_same<T, float>::value) {
        arm_float_to_double(x, x_double, num_samples);
    } else if constexpr (std::is_same<T, q31_t>::value) {
        arm_q31_to_f64(x, x_double, num_samples);
    } else if constexpr (std::is_same<T, double>::value) {
        std::copy(x, x + num_samples, x_double);
    } else {
        static_assert(false, "Unsupported conversion type");
    }
}

template<typename T, typename U>
static void to_native(const U *x, T *x_native, uint32_t num_samples) {
    if constexpr (std::is_same<T, float>::value) {
        if constexpr (std::is_same<U, float>::value) {
            std::copy(x, x + num_samples, x_native);
        } else if constexpr (std::is_same<U, double>::value) {
            arm_float_to_f64(x, x_native, num_samples);
        } else if constexpr (std::is_same<T, q31_t>::value) {
            arm_q31_to_f64(x, x_native, num_samples);
        } else {
            static_assert(false, "Unsupported conversion type");
        }
    } else if constexpr (std::is_same<T, double>::value) {
        if constexpr (std::is_same<U, float>::value) {
            arm_f64_to_float(x, x_native, num_samples);
        } else if constexpr (std::is_same<U, double>::value) {
            std::copy(x, x + num_samples, x_native);
        } else if constexpr (std::is_same<T, q31_t>::value) {
            arm_f64_to_q31(x, x_native, num_samples);
        } else {
            static_assert(false, "Unsupported conversion type");
        }
    } else if constexpr (std::is_same<T, q31_t>::value) {
        if constexpr (std::is_same<U, float>::value) {
            arm_float_to_q31(x, x_native, num_samples);
        } else if constexpr (std::is_same<U, double>::value) {
            arm_f64_to_q31(x, x_native, num_samples);
        } else if constexpr (std::is_same<T, U>::value) {
            std::copy(x, x + num_samples, x_native);
        } else {
            static_assert(false, "Unsupported conversion type");
        }
    } else {
        static_assert(false, "Unsupported conversion type");
    }
}

}
