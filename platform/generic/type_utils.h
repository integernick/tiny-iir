#pragma once

#include <algorithm>
#include <cinttypes>
#include <string>

namespace tiny_iir {

template<typename>
struct always_false : std::false_type {
};

template<typename T>
static void to_double(const T *x, double *x_double, uint32_t num_samples) {
    if constexpr (std::is_same_v<T, float>) {
        std::transform(x, x + num_samples, x_double, [](float sample) {
            return static_cast<double>(sample);
        });
    } else if constexpr (std::is_same_v<T, double>) {
        std::copy(x, x + num_samples, x_double);
    } else {
        static_assert(always_false<T>::value, "Unsupported conversion type");
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
            static_assert(always_false<T>::value, "Unsupported conversion type");
        }
    } else if constexpr (std::is_same_v<T, double>) {
        if constexpr (std::is_same_v<U, float>) {
            std::transform(x, x + num_samples, x_native, [](float sample) {
                return static_cast<double>(sample);
            });
        } else if constexpr (std::is_same_v<U, double>) {
            std::copy(x, x + num_samples, x_native);
        } else {
            static_assert(always_false<T>::value, "Unsupported conversion type");
        }
    } else {
        static_assert(always_false<T>::value, "Unsupported conversion type");
    }
}

template<typename T>
static void multiply(T *x, T *y, T *dst, uint32_t num_samples) {
    if constexpr (std::is_same_v<T, float>) {
        std::transform(x, x + num_samples, y, dst, [](float x, float y) {
            return x * y;
        });
    } else if constexpr (std::is_same_v<T, double>) {
        std::transform(x, x + num_samples, y, dst, [](double x, double y) {
            return x * y;
        });
    } else {
        static_assert(always_false<T>::value, "Unsupported conversion type");
    }
}

template<typename T>
static void scale(const T *src, const T &scale, T *dst, uint32_t num_samples) {
    if constexpr (std::is_same_v<T, float>) {
        std::transform(src, src + num_samples, dst, [&scale](float x) {
            return x * scale;
        });
    } else if constexpr (std::is_same_v<T, double>) {
        std::transform(src, src + num_samples, dst, [&scale](double x) {
            return x * scale;
        });
    } else {
        static_assert(always_false<T>::value, "Unsupported conversion type");
    }
}

} // namespace tiny_iir
