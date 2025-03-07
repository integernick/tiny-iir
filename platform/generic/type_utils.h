#pragma once

#include <algorithm>
#include <string>

namespace tiny_iir {

template<typename T>
static void to_double(const T *x, double *x_double, uint32_t num_samples) {
    if constexpr (std::is_same_v<T, float>) {
        std::transform(x, x + num_samples, x_double, [](float sample) {
            return static_cast<double>(sample);
        });
    } else if constexpr (std::is_same_v<T, double>) {
        std::copy(x, x + num_samples, x_double);
    } else {
        static_assert(false, "Unsupported conversion type");
    }
}

template<typename T, typename U>
static void to_native(const U *x, T *x_native, uint32_t num_samples) {
    if constexpr (std::is_same<T, float>::value) {
        if constexpr (std::is_same_v<U, float>) {
            std::copy(x, x + num_samples, x_native);
        } else if constexpr (std::is_same_v<U, double>) {
            std::transform(x, x + num_samples, x_native, [](double sample) {
                return static_cast<float>(sample);
            });
        } else {
            static_assert(false, "Unsupported conversion type");
        }
    } else if constexpr (std::is_same_v<T, double>) {
        if constexpr (std::is_same_v<U, float>) {
            arm_f64_to_float(x, x_native, num_samples);
        } else if constexpr (std::is_same_v<U, double>) {
            std::copy(x, x + num_samples, x_native);
        } else {
            static_assert(false, "Unsupported conversion type");
        }
    } else {
        static_assert(false, "Unsupported conversion type");
    }
}

}
