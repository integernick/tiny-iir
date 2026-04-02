#pragma once

#include <cstdint>
#include <cstring>
#include <type_traits>

namespace tiny_iir {

template<typename>
struct always_false : std::false_type {
};

template<typename T>
static void to_double(const T *x, double *x_double, uint32_t num_samples) {
    if constexpr (std::is_same_v<T, float>) {
        for (uint32_t i = 0; i < num_samples; ++i) {
            x_double[i] = static_cast<double>(x[i]);
        }
    } else if constexpr (std::is_same_v<T, double>) {
        memcpy(x_double, x, num_samples * sizeof(double));
    } else {
        static_assert(always_false<T>::value, "Unsupported conversion type");
    }
}

template<typename T, typename U>
static void to_native(const U *x, T *x_native, uint32_t num_samples) {
    if constexpr (std::is_same<T, float>::value) {
        if constexpr (std::is_same_v<U, float>) {
            memcpy(x_native, x, num_samples * sizeof(float));
        } else if constexpr (std::is_same_v<U, double>) {
            for (uint32_t i = 0; i < num_samples; ++i) {
                x_native[i] = static_cast<float>(x[i]);
            }
        } else {
            static_assert(always_false<T>::value, "Unsupported conversion type");
        }
    } else if constexpr (std::is_same_v<T, double>) {
        if constexpr (std::is_same_v<U, float>) {
            for (uint32_t i = 0; i < num_samples; ++i) {
                x_native[i] = static_cast<double>(x[i]);
            }
        } else if constexpr (std::is_same_v<U, double>) {
            memcpy(x_native, x, num_samples * sizeof(double));
        } else {
            static_assert(always_false<T>::value, "Unsupported conversion type");
        }
    } else {
        static_assert(always_false<T>::value, "Unsupported conversion type");
    }
}

template<typename T>
static void multiply(T *x, T *y, T *dst, uint32_t num_samples) {
    for (uint32_t i = 0; i < num_samples; ++i) {
        dst[i] = x[i] * y[i];
    }
}

template<typename T>
static void scale(const T *src, const T &scale_val, T *dst, uint32_t num_samples) {
    for (uint32_t i = 0; i < num_samples; ++i) {
        dst[i] = src[i] * scale_val;
    }
}

} // namespace tiny_iir
