#pragma once

#include <cstdint>

namespace tiny_iir {
/**
 * @addtogroup tiny_iir_module
 * @{
 */

template<template<typename> class Cascade, typename T = double>
void biquad_cascade_init(Cascade<T> *cascade, uint8_t num_stages,
                         T *coefficients, T *delay_pipeline);

template<template<typename> class Cascade, typename T = double>
void biquad_cascade_process(Cascade<T> *cascade, T *src, T *dst, uint32_t block_size);

template<template<typename> class Cascade, typename T>
void biquad_cascade_init(Cascade<T> *cascade, uint8_t num_stages,
                         T *coefficients, T *delay_pipeline) {
    if (cascade) {
        cascade->init(num_stages, coefficients, delay_pipeline);
    }
}

template<template<typename> class Cascade, typename T>
void biquad_cascade_process(Cascade<T> *cascade, T *src, T *dst, uint32_t block_size) {
    if (cascade) {
        cascade->process(src, dst, block_size);
    }
}

static constexpr unsigned int COEFFICIENTS_PER_BIQUAD_BLOCK = 5;

/**
 * @brief  Biquad block – (Direct Form I).
 *
 * @details
 * Implements the transfer function
 * \f[
 *     H(z)=\frac{b_0+b_1 z^{-1}+b_2 z^{-2}}
 *                {1-a_1 z^{-1}-a_2 z^{-2}}
 * \f]
 * using the four‑state DF‑I structure:
 * \f[
 * \begin{aligned}
 * y[n] &= b_0\,x[n] + b_1\,x[n-1] + b_2\,x[n-2] \\
 *      &\quad + a_1\,y[n-1] + a_2\,y[n-2] \\[6pt]
 * x[n-2] &\leftarrow x[n-1],\qquad x[n-1] \leftarrow x[n] \\
 * y[n-2] &\leftarrow y[n-1],\qquad y[n-1] \leftarrow y[n]
 * \end{aligned}
 * \f]
 * which requires **four** delay elements (\c x1, \c x2, \c y1, \c y2) per stage.
 *
 * @tparam T  Data type.
 */
template<typename T = double>
class BiquadBlockDF1Cascade {
public:
    void init(uint8_t num_stages, T *coefficients, T *delay_pipeline) {
        _num_stages = num_stages;
        _coefficients = coefficients;
        _delay_pipeline = delay_pipeline;
    }

    void process(T *src, T *dst, uint32_t block_size) {
        for (uint32_t i = 0; i < block_size; ++i) {
            T val = src[i];

            for (uint32_t stage_idx = 0; stage_idx < _num_stages; ++stage_idx) {
                val = process(val, stage_idx);
            }

            dst[i] = val;
        }
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 4;

private:
    T process(T input, uint32_t stage_idx) {
        T *B = _coefficients + COEFFICIENTS_PER_BIQUAD_BLOCK * stage_idx;
        T *A = B + 3;
        T *delay_line_x = _delay_pipeline + BLOCK_DELAY_LINE_SIZE * stage_idx;
        T *delay_line_y = delay_line_x + 2;

        const T output = B[0] * input + B[1] * delay_line_x[0] + B[2] * delay_line_x[1]
                         + A[0] * delay_line_y[0] + A[1] * delay_line_y[1];

        delay_line_x[1] = delay_line_x[0];
        delay_line_x[0] = input;
        delay_line_y[1] = delay_line_y[0];
        delay_line_y[0] = output;

        return output;
    }

    T *_coefficients;
    T *_delay_pipeline;

    uint32_t _num_stages;
};

/**
 * @brief   Biquad block (Direct Form II transposed).
 *
 * @details
 * Implements the transfer function
 * \f[
 *     H(z)=\frac{b_0+b_1 z^{-1}+b_2 z^{-2}}
 *                {1-a_1 z^{-1}-a_2 z^{-2}}
 * \f]
 * using the two‑state DF‑II‑T structure:
 * \f[
 * \begin{aligned}
 * y[n] &= b_0\,x[n] + d_1[n] \\
 * s_1[n+1] &= b_1\,x[n] + a_1\,y[n] + s_2[n] \\
 * s_2[n+1] &= b_2\,x[n] + a_2\,y[n]
 * \end{aligned}
 * \f]
 * which requires only two delay elements \c d1 and \c d2 per stage.
 *
 * @tparam T  Data type.
 */
template<typename T = double>
class BiquadBlockDF2TCascade {
public:
    void init(uint8_t num_stages, T *coefficients, T *delay_pipeline) {
        _num_stages = num_stages;
        _coefficients = coefficients;
        _delay_pipeline = delay_pipeline;
    }

    void process(T *src, T *dst, uint32_t block_size) {
        for (uint32_t i = 0; i < block_size; ++i) {
            T val = src[i];

            for (uint32_t stage_idx = 0; stage_idx < _num_stages; ++stage_idx) {
                val = process(val, stage_idx);
            }

            dst[i] = val;
        }
    }

    static constexpr uint32_t BLOCK_DELAY_LINE_SIZE = 2;

private:
    T process(T input, uint32_t stage_idx) {
        T *B = _coefficients + COEFFICIENTS_PER_BIQUAD_BLOCK * stage_idx;
        T *A = B + 3;
        T *s = _delay_pipeline + BLOCK_DELAY_LINE_SIZE * stage_idx;

        T output = B[0] * input + s[0];
        s[0] = B[1] * input + A[0] * output + s[1];
        s[1] = B[2] * input + A[1] * output;

        return output;
    }

    T *_coefficients;
    T *_delay_pipeline;

    uint32_t _num_stages;
};

/**
 * @}
 */

}