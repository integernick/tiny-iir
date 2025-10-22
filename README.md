# tiny–iir  
small header-only iir filters for realtime / embedded

## what it is
- butterworth / chebyshev i & ii / elliptic  
- low / high / band‑pass / band‑stop  
- cmsis-dsp or generic
- float/double; q15/q31 (cmsis only)

## include

### cmsis‑dsp
```cmake
set(TINY_IIR_CMSIS_CORE_DIR path/to/CMSIS/Include)
add_subdirectory(path/to/tiny-iir)
target_link_libraries(my-project PRIVATE tiny_iir_core)
```

### no cmsis‑dsp
```cmake
set(TINY_IIR_BUILD_WITH_CMSIS OFF CACHE BOOL "Build using CMSIS-DSP")
add_subdirectory(path/to/tiny-iir)
target_link_libraries(my-project PRIVATE tiny_iir_core)
```

### chunk size
`tiny‑iir` uses a small on‑stack scratch buffer for block processing.  
default is **32**. tune via a compile‑time define:
```cmake
# choose bounded on‑stack chunk size (default 32)
target_compile_definitions(my-project PRIVATE TINY_IIR_CHUNK_SIZE=128)
```

## use


### make a filter

**template parameters**
- `order` = analog prototype order (typically 1–20)
- `runtime type` = float, double, q15 or q31
- `pass type`

**arguments**
- frequencies are normalized to nyquist: `1.0 == fs/2`
  - low/high‑pass: single `cutoff` in `(0, 1)`
  - band‑pass / band‑stop: `lowcut`, `highcut` with `0 < lowcut < highcut < 1`
- ripple / attenuation in dB
  - cheby1: passband ripple in dB
  - cheby2: stopband attenuation in dB
  - elliptic: passband ripple and stopband attenuation, both in dB
- (optional) number of crossfade samples to smooth reconfig transitions; `0` disables (default)

```cpp
#include <cheby2/IIRCheby2.h>

// order 6, low‑pass, stopband 60 dB, 100 crossfade samples
tiny_iir::IIRCheby2<6, double, tiny_iir::FilterPassType::LowPass> f{
    0.1, 60.0, 100
};
```

```cpp
#include <elliptic/IIRElliptic.h>

// order 9, band‑pass [0.4, 0.6], 0.1 dB ripple, 60 dB stopband, 1000 crossfade samples
tiny_iir::IIRElliptic<9, float, tiny_iir::FilterPassType::BandPass> g{
    0.4f, 0.6f, 0.1f, 60.0f, 1000
};
```

### process
```cpp
// single sample
double y = f.process(1.0);

// batch — returns last sample
constexpr size_t n = 4;
double x[n]{0.1, 0.2, 0.3, 0.4};
double last = f.process(x, n);
```

### tune
```cpp
// re‑configure cutoff / ripple
f.configure(0.2, 50.0);

// optionally reset internal state (delay lines)
f.reset_state();
```

## design precision
design precision is a template knob. runtime/sample type stays **T**.  
pick **float** for smaller code + faster reconfig; **double** for max stability.

```cpp
// runtime q31, design math in float
tiny_iir::IIRCheby2<6, q31_t, tiny_iir::FilterPassType::LowPass, float> h{0.1f, 60.0f};
```

## cli
generate sos as json for plotting or shipping:
```sh
tiny-iir-designer-cli -t cheby1 -p bsf -o 13 -l 0.2 -h 0.9 --ripple 0.01 --stop 40 -j out.json
```

## plot
simple python tool to preview frequency response:
```sh
python3 tools/plot_sos.py [path/to/coeff_presets.json] [preset_key]
```

examples:
```sh
python3 tools/plot_sos.py tools/demo.json butter_lpf     # butterworth low‑pass
```
![butterworth low‑pass](figures/butter_lpf.png)

```sh
python3 tools/plot_sos.py tools/demo.json cheby1_hpf     # chebyshev i high‑pass
```
![chebyshev i high‑pass](figures/cheby1_hpf.png)

```sh
python3 tools/plot_sos.py tools/demo.json cheby2_bpf     # chebyshev ii band‑pass
```
![chebyshev ii band‑pass](figures/cheby2_bpf.png)

```sh
python3 tools/plot_sos.py tools/demo.json elliptic_bsf   # elliptic band‑stop
```
![elliptic band‑stop](figures/elliptic_bsf.png)
