#include "circuit.hpp"
#include <complex>
#include <string>


std::vector<std::complex<double>> amplitudes(qubit& q) {
    std::complex<double> alpha = std::cos(q.theta() / 2.0);
    std::complex<double> beta  = std::sin(q.theta() / 2.0) * std::exp(std::complex<double>(0.0, q.phi()));
    return {alpha, beta};
}

std::string print_amplitudes(qubit& q) {
    auto amps = amplitudes(q);

    std::string result;

    for (int i = 0; i < 2; ++i) {
        double re = std::round(amps[i].real() * 100.0) / 100.0;
        double im = std::round(amps[i].imag() * 100.0) / 100.0;

        // Convert to string with 2 decimals manually
        int re_int = static_cast<int>(re);
        int re_frac = static_cast<int>(std::round(std::abs(re - re_int) * 100));
        std::string re_str = std::to_string(re_int) + "." + (re_frac < 10 ? "0" : "") + std::to_string(re_frac);

        std::string part = "(" + re_str;

        // Only include imaginary part if non-zero
        if (im != 0.0) {
            std::string sign = (im >= 0 ? " + " : " - ");
            int im_int = static_cast<int>(std::abs(im));
            int im_frac = static_cast<int>(std::round((std::abs(im) - im_int) * 100));
            std::string im_str = std::to_string(im_int) + "." + (im_frac < 10 ? "0" : "") + std::to_string(im_frac);
            part += sign + im_str + "i";
        }

        part += ")";
        part += (i == 0 ? "|0> + " : "|1>");
        result += part;
    }

    return result;
}


/*
Simulate the probability of measuring a qubit in state |1> after applying a circuit of gates.
*/
double simulate_probability(int num_runs, circuit c, qubit& q) {
    int count_one = 0;

    for (int i = 0; i < num_runs; ++i) {
        q.reset();
        c.apply_to(q);
        if (q.measure() == 1) {
            ++count_one;
        }
    }

    return static_cast<double>(count_one) / num_runs;
}
