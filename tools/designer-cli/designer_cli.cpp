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
std::string sos_line(const BiquadCoefficients &c) {
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
    auto invoke = [&](auto N_c) {
        constexpr int N = N_c.value;
        Family<N, double, P> iir(std::forward<Args>(args)...);
        gain_out = iir.get_gain();
        collect<Family<N, double, P>::NUMBER_OF_BIQUAD_BLOCKS>(iir, sos_out);
    };

    switch (order) {
        case 1:
            invoke(std::integral_constant<int, 1>{});
            break;
        case 2:
            invoke(std::integral_constant<int, 2>{});
            break;
        case 3:
            invoke(std::integral_constant<int, 3>{});
            break;
        case 4:
            invoke(std::integral_constant<int, 4>{});
            break;
        case 5:
            invoke(std::integral_constant<int, 5>{});
            break;
        case 6:
            invoke(std::integral_constant<int, 6>{});
            break;
        case 7:
            invoke(std::integral_constant<int, 7>{});
            break;
        case 8:
            invoke(std::integral_constant<int, 8>{});
            break;
        case 9:
            invoke(std::integral_constant<int, 9>{});
            break;
        case 10:
            invoke(std::integral_constant<int, 10>{});
            break;
        case 11:
            invoke(std::integral_constant<int, 11>{});
            break;
        case 12:
            invoke(std::integral_constant<int, 12>{});
            break;
        case 13:
            invoke(std::integral_constant<int, 13>{});
            break;
        case 14:
            invoke(std::integral_constant<int, 14>{});
            break;
        case 15:
            invoke(std::integral_constant<int, 15>{});
            break;
        case 16:
            invoke(std::integral_constant<int, 16>{});
            break;
        case 17:
            invoke(std::integral_constant<int, 17>{});
            break;
        case 18:
            invoke(std::integral_constant<int, 18>{});
            break;
        case 19:
            invoke(std::integral_constant<int, 19>{});
            break;
        case 20:
            invoke(std::integral_constant<int, 20>{});
            break;
        default:
            throw std::invalid_argument("order 1‑20 only");
    }
}

template<int NUM_OF_BIQUAD_BLOCKS, typename T>
void collect(T &f, std::vector<std::string> &v) {
    const auto *s = reinterpret_cast<const BiquadCoefficients *>(f.get_coefficients());

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
            {"lpf", FilterPassType::LOW_PASS},
            {"hpf", FilterPassType::HIGH_PASS},
            {"bpf", FilterPassType::BAND_PASS},
            {"bsf", FilterPassType::BAND_STOP},
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

    if ((pass == FilterPassType::BAND_PASS || pass == FilterPassType::BAND_STOP)) {
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
        if (pass == FilterPassType::LOW_PASS) {
            make_filter<IIRButter, FilterPassType::LOW_PASS>(order, sos, gain, fc1);
        } else if (pass == FilterPassType::HIGH_PASS) {
            make_filter<IIRButter, FilterPassType::HIGH_PASS>(order, sos, gain, fc1);
        } else if (pass == FilterPassType::BAND_PASS) {
            make_filter<IIRButter, FilterPassType::BAND_PASS>(order, sos, gain, fc1, fc2);
        } else {
            make_filter<IIRButter, FilterPassType::BAND_STOP>(order, sos, gain, fc1, fc2);
        }
    } else if (type == "cheby1") {
        if (pass == FilterPassType::LOW_PASS) {
            make_filter<IIRCheby1, FilterPassType::LOW_PASS>(order, sos, gain, fc1, rp);
        } else if (pass == FilterPassType::HIGH_PASS) {
            make_filter<IIRCheby1, FilterPassType::HIGH_PASS>(order, sos, gain, fc1, rp);
        } else if (pass == FilterPassType::BAND_PASS) {
            make_filter<IIRCheby1, FilterPassType::BAND_PASS>(order, sos, gain, fc1, fc2, rp);
        } else {
            make_filter<IIRCheby1, FilterPassType::BAND_STOP>(order, sos, gain, fc1, fc2, rp);
        }
    } else if (type == "cheby2") {
        if (pass == FilterPassType::LOW_PASS) {
            make_filter<IIRCheby2, FilterPassType::LOW_PASS>(order, sos, gain, fc1, rs);
        } else if (pass == FilterPassType::HIGH_PASS) {
            make_filter<IIRCheby2, FilterPassType::HIGH_PASS>(order, sos, gain, fc1, rs);
        } else if (pass == FilterPassType::BAND_PASS) {
            make_filter<IIRCheby2, FilterPassType::BAND_PASS>(order, sos, gain, fc1, fc2, rs);
        } else {
            make_filter<IIRCheby2, FilterPassType::BAND_STOP>(order, sos, gain, fc1, fc2, rs);
        }
    } else if (type == "elliptic") {
        if (pass == FilterPassType::LOW_PASS) {
            make_filter<IIRElliptic, FilterPassType::LOW_PASS>(order, sos, gain, fc1, rp, rs);
        } else if (pass == FilterPassType::HIGH_PASS) {
            make_filter<IIRElliptic, FilterPassType::HIGH_PASS>(order, sos, gain, fc1, rp, rs);
        } else if (pass == FilterPassType::BAND_PASS) {
            make_filter<IIRElliptic, FilterPassType::BAND_PASS>(order, sos, gain, fc1, fc2, rp, rs);
        } else {
            make_filter<IIRElliptic, FilterPassType::BAND_STOP>(order, sos, gain, fc1, fc2, rp, rs);
        }
    } else {
        std::cerr << "Supported filter families:\n\tbutter cheby1 cheby2 elliptic\n";
        return 1;
    }

    const std::string key = type + "_" + pass_str;

    json out;
    json cfg;

    cfg["order"] = order;

    if (pass == FilterPassType::LOW_PASS || pass == FilterPassType::HIGH_PASS) {
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
