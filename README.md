# tiny-iir
Tiny-IIR is a high-performance, floating- and fixed-point Infinite Impulse Response (IIR)
filter library designed for embedded systems and real-time applications using ARM Cortex-M microcontrollers.

## Including the library with CMSISDSP

To use the library with CMSISDSP, simply add the following lines to your `CMakeLists.txt` file:

```cmake
set(TINY_IIR_CMSIS_CORE_DIR path/to/CMSIS/Include)
add_subdirectory(path/to/tiny-iir)
target_link_libraries(my-project PRIVATE tiny_iir_core)
```

## Including the library without CMSISDSP

To include the library in your project, simply add the following lines to your `CMakeLists.txt` file:

```cmake
set(BUILD_WITH_CMSIS OFF CACHE BOOL "Build using CMSIS-DSP")
add_subdirectory(path/to/tiny-iir)
target_link_libraries(my-project PRIVATE tiny_iir_core)
```

To design a filter simply include a corresponding header file and instantiate the filter object.

```cpp
#include <cheby2/IIRCheby2.h>

/* Create an order 6 Chebyshev Type II filter with a passband of 0.1, 
   a stopband of 60dB and 100 crossfade samples */
tiny_iir::IIRCheby2<6, double, tiny_iir::FilterPassType::LOW_PASS> iir_cheby2{
    0.1, 60.0, 100
};
```

```cpp
#include <elliptic/IIRElliptic.h>

/* Create an order 9 elliptic filter with a passband of [0.4, 0.6], 
   passband ripple of 0.1dB, stopband of 60dB (no crossfade) */
tiny_iir::IIRElliptic<9, double, tiny_iir::FilterPassType::BAND_PASS> iir_elliptic{
    0.4, 0.6, 0.1, 60.0
};
```

To process new samples, simply call the `process` method with the new sample as an argument.

```cpp
/* Process a single sample */
double input = 1.0;
const double output = iir_cheby2.process(input);
```

To process a batch of samples, simply call the `process` method with the new sample as an argument.

```cpp
/* Process a batch of samples */
constexpr size_t NUM_SAMPLES = 4;
const double input[NUM_SAMPLES] = {1.0, 2.0, 3.0, 4.0};
double output = iir_cheby2.process(input, NUM_SAMPLES);
```

## Advanced Usage

The filter class provides a number of methods to configure the filter, reset the filter state, and access the filter coefficients.
```cpp
/* Reconfigure the previously created Chebyshev Type II filter 
   to a different cutoff frequency 0.2 and stopband ripple of 50dB */
iir_cheby2.configure(0.2, 50.0);
```

```cpp
/* Reset the filter state */
iir_cheby2.reset_state();
```

## Plotting tool

A simple Python script is provided to plot the frequency response of the filter.
To use it, simply run the following command:

```sh
python3 tools/plot_sos.py [path/to/coeff_presets.json] [preset_key]
```

A number of presets are provided in the `tools/demo.json` file.
For example, to plot the frequency response of the Chebyshev Type I high-pass filter from the `demo.json` file,
run the following command:

```sh
python3 tools/plot_sos.py tools/demo.json cheby1_hpf
```