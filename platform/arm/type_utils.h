#pragma once

#include <dsp/filtering_functions.h>
#include <string>

namespace tiny_iir {

template <typename>
struct always_false : std::false_type {
};

/**
 * @brief   Convert array to double in-place.
 *
 * @note    This function is only available for float and double types.
 * @param x  The input array.
 * @param x_double  The output array.
 * @param num_samples  The number of samples.
 */
template<typename T>
static void to_double(const T *x, double *x_double, uint32_t num_samples) {
    if constexpr (std::is_same_v<T, float>) {
        arm_float_to_f64(x, x_double, num_samples);
    } else if constexpr (std::is_same_v<T, q31_t>) {
        arm_q31_to_f64(x, x_double, num_samples);
    } else if constexpr (std::is_same_v<T, double>) {
        std::copy(x, x + num_samples, x_double);
    } else if constexpr (std::is_same_v<T, q15_t>) {
        arm_q15_to_f64(x, x_double, num_samples);
    } else {
        static_assert(always_false<T>::value, "Unsupported conversion type");
    }
}

/**
 * @brief   Convert array to native type in-place.
 *
 * @note    This function is only available for float and double types.
 * @param x  The input array.
 * @param x_native  The output array of native type.
 * @param num_samples  The number of samples.
 */
template<typename T, typename U>
static void to_native(const U *x, T *x_native, uint32_t num_samples) {
    if constexpr (std::is_same_v<T, float>) {
        if constexpr (std::is_same_v<U, float>) {
            std::copy(x, x + num_samples, x_native);
        } else if constexpr (std::is_same_v<U, double>) {
            arm_f64_to_float(x, x_native, num_samples);
        } else if constexpr (std::is_same_v<U, q31_t>) {
            arm_q31_to_float(x, x_native, num_samples);
        } else if constexpr (std::is_same_v<U, q15_t>) {
            arm_q15_to_float(x, x_native, num_samples);
        } else {
            static_assert(always_false<T>::value, "Unsupported conversion type");
        }
    } else if constexpr (std::is_same_v<T, double>) {
        if constexpr (std::is_same_v<U, float>) {
            arm_float_to_f64(x, x_native, num_samples);
        } else if constexpr (std::is_same_v<U, double>) {
            std::copy(x, x + num_samples, x_native);
        } else if constexpr (std::is_same_v<T, q31_t>) {
            arm_q31_to_f64(x, x_native, num_samples);
        } else if constexpr (std::is_same_v<T, q15_t>) {
            arm_q15_to_f64(x, x_native, num_samples);
        } else {
            static_assert(always_false<T>::value, "Unsupported conversion type");
        }
    } else if constexpr (std::is_same_v<T, q31_t>) {
        if constexpr (std::is_same_v<U, float>) {
            arm_float_to_q31(x, x_native, num_samples);
        } else if constexpr (std::is_same_v<U, double>) {
            arm_f64_to_q31(x, x_native, num_samples);
        } else if constexpr (std::is_same_v<T, U>) {
            std::copy(x, x + num_samples, x_native);
        } else {
            static_assert(always_false<T>::value, "Unsupported conversion type");
        }
    } else if constexpr (std::is_same_v<T, q15_t>) {
        if constexpr (std::is_same_v<U, float>) {
            arm_float_to_q15(x, x_native, num_samples);
        } else if constexpr (std::is_same_v<U, double>) {
            arm_f64_to_q15(x, x_native, num_samples);
        } else if constexpr (std::is_same_v<T, U>) {
            std::copy(x, x + num_samples, x_native);
        } else {
            static_assert(always_false<T>::value, "Unsupported conversion type");
        }
    } else {
        static_assert(always_false<T>::value, "Unsupported conversion type");
    }
}

}
