#include "circuit.hpp"
#include <vector>
#include <string>
#include <unordered_map>
#include <random>
#include <cmath>
#include <algorithm>

using StateProb = std::pair<std::string, double>;
using ProbDist = std::vector<StateProb>;

class quantum_computer {
private:
    circuit& qcircuit;

public:
    quantum_computer(circuit& c) : qcircuit(c) {}

    std::string print_prob_dist(const ProbDist& dist, int precision = 2) {
        std::string result;

        double factor = std::pow(10, precision);

        for (const auto& [bits, prob] : dist) {
            // Convert to percentage and round
            double pct = std::round(prob * 100.0 * factor) / factor;

            // Convert to string with fixed decimal places
            std::string pct_str = std::to_string(pct);

            // Remove trailing zeros after decimal
            if (pct_str.find('.') != std::string::npos) {
                // Erase trailing zeros
                pct_str.erase(pct_str.find_last_not_of('0') + 1, std::string::npos);
                // If ends with '.', remove it too
                if (pct_str.back() == '.') pct_str.pop_back();
            }

            result += "|" + bits + "> " + pct_str + "%\n";
        }

        return result;
    }

    // Deterministic probability distribution from amplitudes
    ProbDist calculate_probabilities() {
        ProbDist dist;
        size_t n = qcircuit.size();
        size_t dim = 1ULL << n;

        for (size_t i = 0; i < dim; ++i) {
            double prob = qcircuit.amplitudes[i] * qcircuit.amplitudes[i]; // |amp|^2
            if (prob < 1e-12) continue;

            // Convert index -> bitstring
            std::string bits;
            bits.reserve(n);
            for (size_t q = 0; q < n; ++q) {
                size_t bit = 1ULL << (n - q - 1);
                bits.push_back((i & bit) ? '1' : '0');
            }

            dist.emplace_back(bits, prob);
        }

        return dist;
    }

    // Monte Carlo simulation using amplitudes
    ProbDist simulate(int num_runs) {
        size_t n = qcircuit.size();
        size_t dim = 1ULL << n;
        ProbDist dist;
        std::unordered_map<std::string, int> counts;

        // Precompute cumulative probabilities for measurement
        std::vector<double> cumulative(dim, 0.0);
        cumulative[0] = qcircuit.amplitudes[0] * qcircuit.amplitudes[0];
        for (size_t i = 1; i < dim; ++i) {
            cumulative[i] = cumulative[i - 1] + qcircuit.amplitudes[i] * qcircuit.amplitudes[i];
        }

        // RNG
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dis(0.0, 1.0);

        for (int run = 0; run < num_runs; ++run) {
            double r = dis(gen);

            // Find which basis state corresponds to r
            auto it = std::lower_bound(cumulative.begin(), cumulative.end(), r);
            size_t idx = std::distance(cumulative.begin(), it);

            // Convert index -> bitstring
            std::string bits;
            bits.reserve(n);
            for (size_t q = 0; q < n; ++q) {
                size_t bit = 1ULL << (n - q - 1);
                bits.push_back((idx & bit) ? '1' : '0');
            }

            counts[bits]++;
        }

        // Convert counts -> probabilities
        for (const auto& pair : counts) {
            dist.emplace_back(pair.first, static_cast<double>(pair.second) / num_runs);
        }

        return dist;
    }
};
