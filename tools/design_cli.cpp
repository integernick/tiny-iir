#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#include <nlohmann/json.hpp>
#include <cxxopts.hpp>

#include <butter/IIRButter.h>
#include <cheby1/IIRCheby1.h>
#include <cheby2/IIRCheby2.h>
#include <elliptic/IIRElliptic.h>

using json = nlohmann::json;
using namespace tiny_iir;

static std::string sos_line(const BiquadCoefficients &c) {
    std::ostringstream ss;
    ss.setf(std::ios::fixed);
    ss << std::setprecision(15)
       << 1.0 << ' '
       << c.b1 / c.b0 << ' '
       << c.b2 / c.b0 << ' '
       << 1.0 << ' '
       << -c.a1 / c.b0 << ' '
       << -c.a2 / c.b0;
    return ss.str();
}

static FilterPassType to_pass(const std::string &s) {
    if (s == "lpf" || s == "low" || s == "lowpass") return FilterPassType::LOW_PASS;
    if (s == "hpf" || s == "high" || s == "highpass") return FilterPassType::HIGH_PASS;
    if (s == "bpf" || s == "band" || s == "bandpass") return FilterPassType::BAND_PASS;
    if (s == "brf" || s == "bsf" || s == "bandstop") return FilterPassType::BAND_STOP;
    throw std::invalid_argument("unknown --pass value");
}

#define ORDERS(ORD, CODE)                              \
    switch (ORD) {                                     \
        case  1: CODE( 1); break;                      \
        case  2: CODE( 2); break;                      \
        case  3: CODE( 3); break;                      \
        case  4: CODE( 4); break;                      \
        case  5: CODE( 5); break;                      \
        case  6: CODE( 6); break;                      \
        case  7: CODE( 7); break;                      \
        case  8: CODE( 8); break;                      \
        case  9: CODE( 9); break;                      \
        case 10: CODE(10); break;                      \
        case 11: CODE(11); break;                      \
        case 12: CODE(12); break;                      \
        case 13: CODE(13); break;                      \
        case 14: CODE(14); break;                      \
        case 15: CODE(15); break;                      \
        case 16: CODE(16); break;                      \
        case 17: CODE(17); break;                      \
        case 18: CODE(18); break;                      \
        case 19: CODE(19); break;                      \
        case 20: CODE(20); break;                      \
        default: throw std::invalid_argument("order 1-20 only"); \
    }

template<int NUM_OF_BIQUAD_BLOCKS, typename T>
static void collect(T &f, std::vector<std::string> &v) {
    const auto *s = reinterpret_cast<const BiquadCoefficients *>(f.get_coefficients());
    for (int i = 0; i < NUM_OF_BIQUAD_BLOCKS; ++i) {
        v.push_back(sos_line(s[i]));
    }
}

int main(int argc, char **argv) {
    cxxopts::Options cli("tiny-iir-designer-cli",
                         "Tiny-IIR command-line filter designer\n\n"
                         "  --type : butterworth  cheby1  cheby2  elliptic\n"
                         "  --pass : lpf | hpf | bpf | bsf\n");

    cli.add_options()
            ("t,type", "filter family", cxxopts::value<std::string>())
            ("p,pass", "pass type", cxxopts::value<std::string>())
            ("o,order", "order (1-20)", cxxopts::value<int>())
            ("l,lowcut", "cutoff / low edge", cxxopts::value<double>())
            ("h,highcut", "high edge", cxxopts::value<double>()->default_value("0"))
            ("r,ripple", "passband ripple dB", cxxopts::value<double>()->default_value("0.1"))
            ("s,stop", "stopband atten dB", cxxopts::value<double>()->default_value("40"))
            ("j,json", "output file", cxxopts::value<std::string>())
            ("help", "help");

    auto args = cli.parse(argc, argv);
    if (argc == 1 || args.count("help")) {
        std::cout << cli.help() << '\n';
        return 0;
    }

    const std::string type = args["type"].as<std::string>();
    const FilterPassType pass = to_pass(args["pass"].as<std::string>());
    const int order = args["order"].as<int>();
    const double fc1 = args["lowcut"].as<double>();
    const double fc2 = args["highcut"].as<double>();
    const double rp = args["ripple"].as<double>();
    const double rs = args["stop"].as<double>();

    if (order < 1 || order > 20) {
        std::cerr << "order 1-20\n";
        return 1;
    }
    if (fc1 <= 0 || fc1 >= 1) {
        std::cerr << "lowcut 0-1\n";
        return 1;
    }
    if ((pass == FilterPassType::BAND_PASS || pass == FilterPassType::BAND_STOP) &&
        (fc2 <= fc1 || fc2 >= 1)) {
        std::cerr << "highcut > lowcut and <1\n";
        return 1;
    }

    std::vector<std::string> sos;
    double gain = 1.0;

    /* Dispatch */
    
#define BUTTER_LPF(N) { IIRButter<N,double> f(fc1); gain=f.get_gain(); collect<IIRButter<N,double>::NUMBER_OF_BIQUAD_BLOCKS>(f,sos); }
#define BUTTER_HPF(N) { IIRButter<N,double,FilterPassType::HIGH_PASS> f(fc1); gain=f.get_gain(); collect<IIRButter<N,double,FilterPassType::HIGH_PASS>::NUMBER_OF_BIQUAD_BLOCKS>(f,sos); }
#define BUTTER_BPF(N) { IIRButter<N,double,FilterPassType::BAND_PASS> f(fc1,fc2); gain=f.get_gain(); collect<IIRButter<N,double,FilterPassType::BAND_PASS>::NUMBER_OF_BIQUAD_BLOCKS>(f,sos); }
#define BUTTER_BSF(N) { IIRButter<N,double,FilterPassType::BAND_STOP> f(fc1,fc2); gain=f.get_gain(); collect<IIRButter<N,double,FilterPassType::BAND_STOP>::NUMBER_OF_BIQUAD_BLOCKS>(f,sos); }

#define CHEBY1_LPF(N) { IIRCheby1<N,double> f(fc1,rp); gain=f.get_gain(); collect<IIRButter<N,double>::NUMBER_OF_BIQUAD_BLOCKS>(f,sos); }
#define CHEBY1_HPF(N) { IIRCheby1<N,double,FilterPassType::HIGH_PASS> f(fc1,rp); gain=f.get_gain(); collect<IIRButter<N,double,FilterPassType::HIGH_PASS>::NUMBER_OF_BIQUAD_BLOCKS>(f,sos); }
#define CHEBY1_BPF(N) { IIRCheby1<N,double,FilterPassType::BAND_PASS> f(fc1,fc2,rp); gain=f.get_gain(); collect<IIRButter<N,double,FilterPassType::BAND_PASS>::NUMBER_OF_BIQUAD_BLOCKS>(f,sos); }
#define CHEBY1_BSF(N) { IIRCheby1<N,double,FilterPassType::BAND_STOP> f(fc1,fc2,rp); gain=f.get_gain(); collect<IIRButter<N,double,FilterPassType::BAND_STOP>::NUMBER_OF_BIQUAD_BLOCKS>(f,sos); }

#define CHEBY2_LPF(N) { IIRCheby2<N,double> f(fc1,rs); gain=f.get_gain(); collect<IIRButter<N,double>::NUMBER_OF_BIQUAD_BLOCKS>(f,sos); }
#define CHEBY2_HPF(N) { IIRCheby2<N,double,FilterPassType::HIGH_PASS> f(fc1,rs); gain=f.get_gain(); collect<IIRButter<N,double,FilterPassType::HIGH_PASS>::NUMBER_OF_BIQUAD_BLOCKS>(f,sos); }
#define CHEBY2_BPF(N) { IIRCheby2<N,double,FilterPassType::BAND_PASS> f(fc1,fc2,rs); gain=f.get_gain(); collect<IIRButter<N,double,FilterPassType::BAND_PASS>::NUMBER_OF_BIQUAD_BLOCKS>(f,sos); }
#define CHEBY2_BSF(N) { IIRCheby2<N,double,FilterPassType::BAND_STOP> f(fc1,fc2,rs); gain=f.get_gain(); collect<IIRButter<N,double,FilterPassType::BAND_STOP>::NUMBER_OF_BIQUAD_BLOCKS>(f,sos); }

#define ELLIP_LPF(N)  { IIRElliptic<N,double> f(fc1,rp,rs); gain=f.get_gain(); collect<IIRButter<N,double>::NUMBER_OF_BIQUAD_BLOCKS>(f,sos); }
#define ELLIP_HPF(N)  { IIRElliptic<N,double,FilterPassType::HIGH_PASS> f(fc1,rp,rs); gain=f.get_gain(); collect<IIRButter<N,double,FilterPassType::HIGH_PASS>::NUMBER_OF_BIQUAD_BLOCKS>(f,sos); }
#define ELLIP_BPF(N)  { IIRElliptic<N,double,FilterPassType::BAND_PASS> f(fc1,fc2,rp,rs); gain=f.get_gain(); collect<IIRButter<N,double,FilterPassType::BAND_PASS>::NUMBER_OF_BIQUAD_BLOCKS>(f,sos); }
#define ELLIP_BSF(N)  { IIRElliptic<N,double,FilterPassType::BAND_STOP> f(fc1,fc2,rp,rs); gain=f.get_gain(); collect<IIRButter<N,double,FilterPassType::BAND_STOP>::NUMBER_OF_BIQUAD_BLOCKS>(f,sos); }

    if (type == "butterworth" || type == "butter") {
        if (pass == FilterPassType::LOW_PASS) ORDERS(order, BUTTER_LPF)
        else if (pass == FilterPassType::HIGH_PASS) ORDERS(order, BUTTER_HPF)
        else if (pass == FilterPassType::BAND_PASS) ORDERS(order, BUTTER_BPF)
        else
            ORDERS(order, BUTTER_BSF)
    } else if (type == "cheby1") {
        if (pass == FilterPassType::LOW_PASS) ORDERS(order, CHEBY1_LPF)
        else if (pass == FilterPassType::HIGH_PASS) ORDERS(order, CHEBY1_HPF)
        else if (pass == FilterPassType::BAND_PASS) ORDERS(order, CHEBY1_BPF)
        else
            ORDERS(order, CHEBY1_BSF)
    } else if (type == "cheby2") {
        if (pass == FilterPassType::LOW_PASS) ORDERS(order, CHEBY2_LPF)
        else if (pass == FilterPassType::HIGH_PASS) ORDERS(order, CHEBY2_HPF)
        else if (pass == FilterPassType::BAND_PASS) ORDERS(order, CHEBY2_BPF)
        else
            ORDERS(order, CHEBY2_BSF)
    } else if (type == "elliptic") {
        if (pass == FilterPassType::LOW_PASS) ORDERS(order, ELLIP_LPF)
        else if (pass == FilterPassType::HIGH_PASS) ORDERS(order, ELLIP_HPF)
        else if (pass == FilterPassType::BAND_PASS) ORDERS(order, ELLIP_BPF)
        else
            ORDERS(order, ELLIP_BSF)
    } else {
        std::cerr << "unknown --type\n";
        return 1;
    }

#undef BUTTER_LPF
#undef BUTTER_HPF
#undef BUTTER_BPF
#undef BUTTER_BSF
#undef CHEBY1_LPF
#undef CHEBY1_HPF
#undef CHEBY1_BPF
#undef CHEBY1_BSF
#undef CHEBY2_LPF
#undef CHEBY2_HPF
#undef CHEBY2_BPF
#undef CHEBY2_BSF
#undef ELLIP_LPF
#undef ELLIP_HPF
#undef ELLIP_BPF
#undef ELLIP_BSF

    const std::string tshort = (type == "butterworth" ? "butter" : type);
    const std::string pshort = (pass == FilterPassType::LOW_PASS ? "lpf" :
                                pass == FilterPassType::HIGH_PASS ? "hpf" :
                                pass == FilterPassType::BAND_PASS ? "bpf" : "bsf");
    const std::string key = tshort + "_" + pshort;

    json out;
    json cfg;

    cfg["order"] = order;
    if (pass == FilterPassType::LOW_PASS || pass == FilterPassType::HIGH_PASS)
        cfg["freq_cutoff"] = fc1;
    else {
        cfg["freq_lowcut"] = fc1;
        cfg["freq_highcut"] = fc2;
    }

    if (type == "cheby1" || type == "elliptic") cfg["passband_ripple_dB"] = rp;
    if (type == "cheby2" || type == "elliptic") cfg["stopband_ripple_dB"] = rs;

    cfg["gain"] = gain;
    cfg["sos"] = sos;
    out[key] = cfg;

    if (args.count("json"))
        std::ofstream(args["json"].as<std::string>()) << out.dump(2) << '\n';
    else
        std::cout << out.dump(2) << '\n';

    return 0;
}
