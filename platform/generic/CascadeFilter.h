#pragma once

#include "BiquadBlock.h"

#include <array>

#define CASCADE_FILTER_DEBUG    0
#if CASCADE_FILTER_DEBUG > 0
#include <iostream>
#include <iomanip>
#endif

namespace tiny_iir {
template<unsigned int ORDER, typename T = double,
        class BiquadBlockDirectForm = BiquadBlockDF1<T>>
class CascadeFilter {
public:
    static constexpr unsigned int NUMBER_OF_BIQUAD_BLOCKS = (ORDER + 1) / 2;
    static constexpr unsigned int COEFFICIENTS_PER_BIQUAD_BLOCK = 5;
    static constexpr unsigned int NUMBER_OF_COEFFICIENTS
            = NUMBER_OF_BIQUAD_BLOCKS * COEFFICIENTS_PER_BIQUAD_BLOCK;

    [[nodiscard]] const T *get_coefficients() const {
        return _coefficients;
    }

#if CASCADE_FILTER_DEBUG > 0
    void print_coefficients() const {
    double gain_double;
    BiquadCascade<T>::to_double(&_gain, &gain_double, 1);
    std::cout << std::setprecision(15) << "gain: " << gain_double << std::endl;
    for (int i = 0; i < NUMBER_OF_BIQUAD_BLOCKS; ++i) {
        double biquad_coefficients_double[COEFFICIENTS_PER_BIQUAD_BLOCK];
        BiquadCascade<T>::to_double(&_coefficients[COEFFICIENTS_PER_BIQUAD_BLOCK * i],
                                    biquad_coefficients_double, COEFFICIENTS_PER_BIQUAD_BLOCK);
        std::cout << std::setprecision(15)
                  << biquad_coefficients_double[0]
                  << ", " << biquad_coefficients_double[1]
                  << ", " << biquad_coefficients_double[2]
                  << ", 1"
                  << ", " << biquad_coefficients_double[3]
                  << ", " << biquad_coefficients_double[4] << ";"
                  << std::endl;
    }
}
#endif

    bool push_biquad_coefficients(double b0, double b1, double b2, double a1, double a2) {
        if (_num_biquad_blocks_set >= NUMBER_OF_BIQUAD_BLOCKS) {
            return false;
        }
        T *current_coefficients_block = &_coefficients[_num_biquad_blocks_set * COEFFICIENTS_PER_BIQUAD_BLOCK];
        current_coefficients_block[0] = static_cast<T>(b0);
        current_coefficients_block[1] = static_cast<T>(b1);
        current_coefficients_block[2] = static_cast<T>(b2);
        current_coefficients_block[3] = static_cast<T>(-a1);
        current_coefficients_block[4] = static_cast<T>(-a2);
        _num_biquad_blocks_set++;
        return true;
    }

    void init_biquad_cascades() {
        for (unsigned int i = 0; i < NUMBER_OF_BIQUAD_BLOCKS; ++i) {
            _biquad_blocks[i].set_coefficients(_coefficients + i * COEFFICIENTS_PER_BIQUAD_BLOCK);
        }
    }

    void reset_state() {
        for (auto &block: _biquad_blocks) {
            block.reset();
        }
    }

    void reset() {
        _num_biquad_blocks_set = 0;
        memset(_coefficients, 0, sizeof(_coefficients));
        reset_state();
    }

    [[nodiscard]] T get_gain() const {
        return _gain;
    }

    [[nodiscard]] int get_number_of_blocks() const {
        return _num_biquad_blocks_set;
    }

    void set_gain(double gain) {
        if constexpr (std::is_same<T, double>::value) {
            _gain = gain;
        } else {
            to_native(&gain, &_gain, 1);
        }
    }

    T process(T x) {
        x *= _gain;
        for (auto &biquad_block: _biquad_blocks) {
            x = biquad_block.process(x);
        }
        return x;
    }

    void process(const T *x, T *out, uint32_t num_samples) {
        if (num_samples == 0) {
            return;
        }
        for (uint32_t i = 0; i < num_samples; ++i) {
            out[i] = process(x[i]);
        }
    }

    [[nodiscard]] T process(const T *x, uint32_t num_samples) {
        if (num_samples == 0) {
            return 0;
        }
        T out;
        for (uint32_t i = 0; i < num_samples; ++i) {
            out = process(x[i]);
        }
        return out;
    }

private:
    T _gain;
    T _coefficients[NUMBER_OF_COEFFICIENTS];
    BiquadBlockDirectForm _biquad_blocks[NUMBER_OF_BIQUAD_BLOCKS];

    uint32_t _num_biquad_blocks_set = 0;
};
}
