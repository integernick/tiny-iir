#pragma once

#include <common/common_utils.h>

#include <dsp/filtering_functions.h>
#include <cmath>

namespace tiny_iir {

template<typename>
struct always_false : std::false_type {
};

/**
 * @brief   Convert array to double in-place.
 *
 * @note    This function is only available for float and double types.
 * @param src  The input array.
 * @param dst  The output array.
 * @param num_samples  The number of samples.
 */
template<typename T>
static void to_double(const T *src, double *dst, uint32_t num_samples) {
    if constexpr (std::is_same_v<T, float>) {
        arm_float_to_f64(src, dst, num_samples);
    } else if constexpr (std::is_same_v<T, q31_t>) {
        arm_q31_to_f64(src, dst, num_samples);
    } else if constexpr (std::is_same_v<T, double>) {
        memcpy(dst, src, num_samples * sizeof(double));
    } else if constexpr (std::is_same_v<T, q15_t>) {
        arm_q15_to_f64(src, dst, num_samples);
    } else {
        static_assert(always_false<T>::value, "Unsupported conversion type");
    }
}

/**
 * @brief   Convert array to native type in-place.
 *
 * @note    This function is only available for float and double types.
 * @param src  The input array.
 * @param dst  The output array of native type.
 * @param num_samples  The number of samples.
 */
template<typename T, typename U>
static void to_native(const U *src, T *dst, uint32_t num_samples) {
    if constexpr (std::is_same_v<T, float>) {
        if constexpr (std::is_same_v<U, float>) {
            memcpy(src, dst, num_samples * sizeof(float));
        } else if constexpr (std::is_same_v<U, double>) {
            arm_f64_to_float(src, dst, num_samples);
        } else if constexpr (std::is_same_v<U, q31_t>) {
            arm_q31_to_float(src, dst, num_samples);
        } else if constexpr (std::is_same_v<U, q15_t>) {
            arm_q15_to_float(src, dst, num_samples);
        } else {
            static_assert(always_false<T>::value, "Unsupported conversion type");
        }
    } else if constexpr (std::is_same_v<T, double>) {
        if constexpr (std::is_same_v<U, float>) {
            arm_float_to_f64(src, dst, num_samples);
        } else if constexpr (std::is_same_v<U, double>) {
            memcpy(dst, src, num_samples * sizeof(double));
        } else if constexpr (std::is_same_v<T, q31_t>) {
            arm_q31_to_f64(src, dst, num_samples);
        } else if constexpr (std::is_same_v<T, q15_t>) {
            arm_q15_to_f64(src, dst, num_samples);
        } else {
            static_assert(always_false<T>::value, "Unsupported conversion type");
        }
    } else if constexpr (std::is_same_v<T, q31_t>) {
        if constexpr (std::is_same_v<U, float>) {
            arm_float_to_q31(src, dst, num_samples);
        } else if constexpr (std::is_same_v<U, double>) {
            arm_f64_to_q31(src, dst, num_samples);
        } else if constexpr (std::is_same_v<T, U>) {
            memcpy(dst, src, num_samples * sizeof(q31_t));
        } else {
            static_assert(always_false<T>::value, "Unsupported conversion type");
        }
    } else if constexpr (std::is_same_v<T, q15_t>) {
        if constexpr (std::is_same_v<U, float>) {
            arm_float_to_q15(src, dst, num_samples);
        } else if constexpr (std::is_same_v<U, double>) {
            arm_f64_to_q15(src, dst, num_samples);
        } else if constexpr (std::is_same_v<T, U>) {
            memcpy(dst, src, num_samples * sizeof(q15_t));
        } else {
            static_assert(always_false<T>::value, "Unsupported conversion type");
        }
    } else {
        static_assert(always_false<T>::value, "Unsupported conversion type");
    }
}

/**
 * @brief   Multiply two arrays.
 *
 * @param src_a  The first array.
 * @param src_b  The second array.
 * @param dst  The destination array.
 * @param num_samples  The number of samples.
 */
template<typename T>
static void multiply(T *src_a, T *src_b, T *dst, uint32_t num_samples) {
    if constexpr (std::is_same_v<T, float>) {
        arm_mult_f32(src_a, src_b, dst, num_samples);
    } else if constexpr (std::is_same_v<T, double>) {
        arm_mult_f64(src_a, src_b, dst, num_samples);
    } else if constexpr (std::is_same_v<T, q31_t>) {
        arm_mult_q31(src_a, src_b, dst, num_samples);
    } else if constexpr (std::is_same_v<T, q15_t>) {
        return arm_mult_q15(src_a, src_b, dst, num_samples);
    } else {
        static_assert(always_false<T>::value, "Unsupported conversion type");
    }
}

template<typename T>
static void scale(const T *src, const T &scale, T *dst, uint32_t num_samples) {
    if constexpr (std::is_same_v<T, float>) {
        arm_scale_f32(src, scale, dst, num_samples);
    } else if constexpr (std::is_same_v<T, double>) {
        arm_scale_f64(src, scale, dst, num_samples);
    } else if constexpr (std::is_same_v<T, q31_t>) {
        arm_scale_q31(src, scale, 0, dst, num_samples);
    } else if constexpr (std::is_same_v<T, q15_t>) {
        arm_scale_q15(src, scale, 0, dst, num_samples);
    } else {
        static_assert(always_false<T>::value, "Unsupported conversion type");
    }
}

} // namespace tiny_iir