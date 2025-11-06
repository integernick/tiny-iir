#include <butter/IIRButter.h>
#include <cheby1/IIRCheby1.h>
#include <cheby2/IIRCheby2.h>
#include <elliptic/IIRElliptic.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#include <nlohmann/json.hpp>
#include <cxxopts.hpp>

using json = nlohmann::json;
using namespace tiny_iir;

namespace {
std::string sos_line(const BiquadCoefficients<double> &c) {
    std::ostringstream ss;
    ss.setf(std::ios::fixed);
    ss << std::setprecision(15)
       << 1.0 << ' '
       << c.b1 << ' '
       << c.b2 << ' '
       << 1.0 << ' '
       << -c.a1 << ' '
       << -c.a2;
    return ss.str();
}

template<int NUM_OF_BIQUAD_BLOCKS, typename T>
static void collect(T &f, std::vector<std::string> &v);

template<template<int, typename, FilterPassType> class Family,
        FilterPassType P,
        typename... Args>
void make_filter(int order,
                 std::vector<std::string> &sos_out,
                 double &gain_out,
                 Args &&... args) {
    constexpr int MAX_ORDER = 20;

    if (order < 1 || order > MAX_ORDER)
        throw std::invalid_argument("order 1-20 only");

    // Instantiate Family<ORDER, double, P>, grab gain, collect SOS
    auto instantiate_and_collect = [&](auto order_tag) {
        constexpr int ORDER = decltype(order_tag)::value;
        Family<ORDER, double, P> iir(std::forward<Args>(args)...);
        gain_out = iir.get_gain();
        collect<Family<ORDER, double, P>::NUMBER_OF_BIQUAD_BLOCKS>(iir, sos_out);
    };

    bool matched = false;
    // Iterate ORDINAL (0..MAX_ORDER-1); run only when (order == ORDINAL+1)
    [&]<std::size_t... ordinal_idx>(std::index_sequence<ordinal_idx...>) {
        ((order == int(ordinal_idx + 1)
          ? (instantiate_and_collect(std::integral_constant<int, int(ordinal_idx + 1)>{}),
                        matched = true, 0)
          : 0), ... );
    }(std::make_index_sequence<MAX_ORDER>{});

    (void) matched;
}

template<int NUM_OF_BIQUAD_BLOCKS, typename T>
void collect(T &f, std::vector<std::string> &v) {
    const auto *s = reinterpret_cast<const BiquadCoefficients<double> *>(f.get_coefficients());

    for (int i = 0; i < NUM_OF_BIQUAD_BLOCKS; ++i) {
        v.push_back(sos_line(s[i]));
    }
}
} // namespace

int main(int argc, char **argv) {
    cxxopts::Options cli("tiny-iir-designer-cli",
                         "Tiny-IIR command-line filter designer\n\n"
                         "  --type : butter | cheby1 | cheby2 | elliptic\n"
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

    const std::unordered_map<std::string, FilterPassType> str2pass_map = {
            {"lpf", FilterPassType::LowPass},
            {"hpf", FilterPassType::HighPass},
            {"bpf", FilterPassType::BandPass},
            {"bsf", FilterPassType::BandStop},
    };

    const std::string type = args["type"].as<std::string>();
    const std::string pass_str = args["pass"].as<std::string>();

    if (str2pass_map.find(pass_str) == str2pass_map.end()) {
        std::cerr << "Supported filter pass types:\n";
        for (const auto &[key, value]: str2pass_map) {
            std::cerr << "\t" << key << "\n";
        }
        std::cerr << "\n";
        return 1;
    }

    const FilterPassType pass = str2pass_map.at(pass_str);
    const int order = args["order"].as<int>();
    const double fc1 = args["lowcut"].as<double>();
    const double fc2 = args["highcut"].as<double>();
    const double rp = args["ripple"].as<double>();
    const double rs = args["stop"].as<double>();

    if (order < 1 || order > 20) {
        std::cerr << "Only filter orders 1-20 are supported\n";
        return 1;
    }

    if (fc1 <= 0.0 || fc1 >= 1.0) {
        std::cerr << "Lowcut frequency out of range (0.0, 1.0)\n";
        return 1;
    }

    if ((pass == FilterPassType::BandPass || pass == FilterPassType::BandStop)) {
        if (fc1 >= fc2) {
            std::cerr << "Lowcut frequency >= highcut frequency\n";
            return 1;
        }

        if (fc2 >= 1.0) {
            std::cerr << "Highcut frequency out of range (0.0, 1.0)\n";
            return 1;
        }
    }

    std::vector<std::string> sos;
    double gain = 1.0;

    if (type == "butter") {
        if (pass == FilterPassType::LowPass) {
            make_filter<IIRButter, FilterPassType::LowPass>(order, sos, gain, fc1);
        } else if (pass == FilterPassType::HighPass) {
            make_filter<IIRButter, FilterPassType::HighPass>(order, sos, gain, fc1);
        } else if (pass == FilterPassType::BandPass) {
            make_filter<IIRButter, FilterPassType::BandPass>(order, sos, gain, fc1, fc2);
        } else {
            make_filter<IIRButter, FilterPassType::BandStop>(order, sos, gain, fc1, fc2);
        }
    } else if (type == "cheby1") {
        if (pass == FilterPassType::LowPass) {
            make_filter<IIRCheby1, FilterPassType::LowPass>(order, sos, gain, fc1, rp);
        } else if (pass == FilterPassType::HighPass) {
            make_filter<IIRCheby1, FilterPassType::HighPass>(order, sos, gain, fc1, rp);
        } else if (pass == FilterPassType::BandPass) {
            make_filter<IIRCheby1, FilterPassType::BandPass>(order, sos, gain, fc1, fc2, rp);
        } else {
            make_filter<IIRCheby1, FilterPassType::BandStop>(order, sos, gain, fc1, fc2, rp);
        }
    } else if (type == "cheby2") {
        if (pass == FilterPassType::LowPass) {
            make_filter<IIRCheby2, FilterPassType::LowPass>(order, sos, gain, fc1, rs);
        } else if (pass == FilterPassType::HighPass) {
            make_filter<IIRCheby2, FilterPassType::HighPass>(order, sos, gain, fc1, rs);
        } else if (pass == FilterPassType::BandPass) {
            make_filter<IIRCheby2, FilterPassType::BandPass>(order, sos, gain, fc1, fc2, rs);
        } else {
            make_filter<IIRCheby2, FilterPassType::BandStop>(order, sos, gain, fc1, fc2, rs);
        }
    } else if (type == "elliptic") {
        if (pass == FilterPassType::LowPass) {
            make_filter<IIRElliptic, FilterPassType::LowPass>(order, sos, gain, fc1, rp, rs);
        } else if (pass == FilterPassType::HighPass) {
            make_filter<IIRElliptic, FilterPassType::HighPass>(order, sos, gain, fc1, rp, rs);
        } else if (pass == FilterPassType::BandPass) {
            make_filter<IIRElliptic, FilterPassType::BandPass>(order, sos, gain, fc1, fc2, rp, rs);
        } else {
            make_filter<IIRElliptic, FilterPassType::BandStop>(order, sos, gain, fc1, fc2, rp, rs);
        }
    } else {
        std::cerr << "Supported filter families:\n\tbutter cheby1 cheby2 elliptic\n";
        return 1;
    }

    const std::string key = type + "_" + pass_str;

    json out;
    json cfg;

    cfg["order"] = order;

    if (pass == FilterPassType::LowPass || pass == FilterPassType::HighPass) {
        cfg["freq_cutoff"] = fc1;
    } else {
        cfg["freq_lowcut"] = fc1;
        cfg["freq_highcut"] = fc2;
    }

    if (type == "cheby1" || type == "elliptic") {
        cfg["passband_ripple_dB"] = rp;
    }

    if (type == "cheby2" || type == "elliptic") {
        cfg["stopband_ripple_dB"] = rs;
    }

    cfg["gain"] = gain;
    cfg["sos"] = sos;
    out[key] = cfg;

    if (args.count("json")) {
        std::ofstream(args["json"].as<std::string>()) << out.dump(2) << '\n';
    } else {
        std::cout << out.dump(2) << '\n';
    }

    return 0;
}
