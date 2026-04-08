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

    static void sort_prob_dist(ProbDist& dist) {
        std::sort(dist.begin(), dist.end(), [](const StateProb& a, const StateProb& b) {
            return a.first < b.first;
        });
    }

public:
    quantum_computer(circuit& c) : qcircuit(c) {}

    std::string print_circuit() const {
        return qcircuit.print_circuit();
    }

    std::string print_amplitudes(int precision = 2) const {
        return qcircuit.print_amplitudes(precision);
    }

    std::string print_prob_dist(const ProbDist& dist, int precision = 2) {
        std::string result;
        double factor = std::pow(10, precision);
        for (const auto& [bits, prob] : dist) {
            double pct = std::round(prob * 100.0 * factor) / factor;
            std::string pct_str = std::to_string(pct);
            if (pct_str.find('.') != std::string::npos) {
                pct_str.erase(pct_str.find_last_not_of('0') + 1, std::string::npos);
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
            // std::norm gives |a + bi|^2 = a^2 + b^2
            double prob = std::norm(qcircuit.amplitudes[i]);
            if (prob < 1e-12) continue;

            std::string bits;
            bits.reserve(n);
            for (size_t q = 0; q < n; ++q) {
                size_t bit = 1ULL << (n - q - 1);
                bits.push_back((i & bit) ? '1' : '0');
            }
            dist.emplace_back(bits, prob);
        }

        sort_prob_dist(dist);
        return dist;
    }

    // Monte Carlo simulation using amplitudes
    ProbDist simulate(int num_runs) {
        size_t n = qcircuit.size();
        size_t dim = 1ULL << n;
        ProbDist dist;
        std::unordered_map<std::string, int> counts;

        // Precompute cumulative probabilities using std::norm for |amp|^2
        std::vector<double> cumulative(dim, 0.0);
        cumulative[0] = std::norm(qcircuit.amplitudes[0]);
        for (size_t i = 1; i < dim; ++i) {
            cumulative[i] = cumulative[i - 1] + std::norm(qcircuit.amplitudes[i]);
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dis(0.0, 1.0);

        for (int run = 0; run < num_runs; ++run) {
            double r = dis(gen);
            auto it = std::lower_bound(cumulative.begin(), cumulative.end(), r);
            size_t idx = std::distance(cumulative.begin(), it);

            std::string bits;
            bits.reserve(n);
            for (size_t q = 0; q < n; ++q) {
                size_t bit = 1ULL << (n - q - 1);
                bits.push_back((idx & bit) ? '1' : '0');
            }
            counts[bits]++;
        }

        for (const auto& pair : counts) {
            dist.emplace_back(pair.first, static_cast<double>(pair.second) / num_runs);
        }

        sort_prob_dist(dist);
        return dist;
    }
};