#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <climits>

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using cx = std::complex<double>;

enum class gate {
    None,
    I,          // identity
    H,          // Hadamard
    X,          // Pauli-X  (bit flip)
    Z,          // Pauli-Z  (phase flip)
    Rk,         // Phase rotation R_k = diag(1, e^(2*pi*i/2^k))
    Control,    // control dot for CNOT
    Target,     // target + for CNOT
    SubCircuit  // marker for applied sub-circuit
};

// Forward declaration
class circuit;

class qubit {
    private:
        std::vector<gate> gates;
        // Track sub-circuit applications: map column_index -> subcircuit_name
        std::map<size_t, std::string> subCircuitNames;
        // Track Rk values: map column_index -> k
        std::map<size_t, int> rkValues;
        
    public:
        int measured : 1;
        int value    : 1;
        // Index of the qubit in the circuit
        size_t index = 0;
        // Each qubit tracks how many gates have been applied to it
        size_t gateCount = 0;
        // Qubit name
        std::string name;
        
        // Constructors
        qubit(const std::string& n = "") : measured(0), value(0), name(n) {}
        
        void addGate(gate g, const std::string& subCircuitName = "", int rk = 0) {
            gates.push_back(g);
            if (g == gate::SubCircuit && !subCircuitName.empty()) {
                subCircuitNames[gateCount] = subCircuitName;
            }
            if (g == gate::Rk) {
                rkValues[gateCount] = rk;
            }
            gateCount++;
        }
        
        gate getGate(size_t index) const {
            if (index < gates.size()) {
                return gates[index];
            }
            return gate::None;
        }
        
        int getRkValue(size_t col) const {
            auto it = rkValues.find(col);
            return (it != rkValues.end()) ? it->second : 0;
        }
        
        const std::map<size_t, std::string>& getSubCircuitNames() const {
            return subCircuitNames;
        }
};

class circuit {
    private:
        std::vector<qubit*> qubits;
        std::string circuit_name;

        void validate_qubit_in_circuit(const qubit& q) const {
            for (auto* qp : qubits) {
                if (qp == &q) return;
            }
            throw std::runtime_error(
                "Error: Qubit \"" + q.name + "\" is not part of circuit \"" + circuit_name + "\"."
            );
        }

        void align_all_qubits() {
            size_t max_count = 0;
            for (auto* q : qubits)
                max_count = std::max(max_count, q->gateCount);

            for (auto* q : qubits) {
                while (q->gateCount < max_count)
                    q->addGate(gate::I);
            }
        }

        // --- Private Math Helpers ---

        void math_H(size_t qubit_idx) {
            size_t n = qubits.size();
            size_t bit = 1ULL << (n - qubit_idx - 1);
            std::vector<cx> new_amplitudes(amplitudes.size());
            const double s = M_SQRT1_2;

            for (size_t i = 0; i < amplitudes.size(); ++i) {
                if ((i & bit) == 0) {
                    size_t partner = i | bit;
                    new_amplitudes[i]       = s * (amplitudes[i] + amplitudes[partner]);
                    new_amplitudes[partner] = s * (amplitudes[i] - amplitudes[partner]);
                }
            }
            amplitudes = std::move(new_amplitudes);
        }

        void math_X(size_t qubit_idx) {
            size_t n = qubits.size();
            size_t bit = 1ULL << (n - qubit_idx - 1);
            std::vector<cx> new_amplitudes = amplitudes;

            for (size_t i = 0; i < amplitudes.size(); ++i) {
                size_t flipped = i ^ bit;
                new_amplitudes[flipped] = amplitudes[i];
            }
            amplitudes = std::move(new_amplitudes);
        }

        void math_Z(size_t qubit_idx) {
            size_t n = qubits.size();
            size_t bit = 1ULL << (n - qubit_idx - 1);

            for (size_t i = 0; i < amplitudes.size(); ++i) {
                if (i & bit)
                    amplitudes[i] = -amplitudes[i];
            }
        }

        // R_k = diag(1, e^(2*pi*i / 2^k))
        // Only the |1> component of the target qubit picks up the phase.
        void math_Rk(size_t qubit_idx, int k) {
            size_t n = qubits.size();
            size_t bit = 1ULL << (n - qubit_idx - 1);
            double angle = 2.0 * M_PI / (1ULL << k);   // 2π / 2^k
            cx phase(std::cos(angle), std::sin(angle));

            for (size_t i = 0; i < amplitudes.size(); ++i) {
                if (i & bit)
                    amplitudes[i] *= phase;
            }
        }

        void math_CNOT(size_t ctrl_idx, size_t tgt_idx) {
            size_t n = qubits.size();
            size_t control_bit = 1ULL << (n - ctrl_idx - 1);
            size_t target_bit  = 1ULL << (n - tgt_idx  - 1);
            std::vector<cx> new_amplitudes = amplitudes;

            for (size_t i = 0; i < amplitudes.size(); ++i) {
                if (i & control_bit) {
                    size_t flipped = i ^ target_bit;
                    new_amplitudes[flipped] = amplitudes[i];
                } else {
                    new_amplitudes[i] = amplitudes[i];
                }
            }
            amplitudes = std::move(new_amplitudes);
        }

    public:
        std::vector<cx> amplitudes;

        circuit(const std::vector<qubit*>& q, const std::string& name = "")
            : qubits(q), circuit_name(name) {
            amplitudes.resize(1 << qubits.size(), cx(0.0, 0.0));
            amplitudes[0] = cx(1.0, 0.0);
            for (size_t i = 0; i < qubits.size(); ++i) {
                qubits[i]->index = i;
            }
        }

        size_t size() const { return qubits.size(); }

        qubit& operator[](size_t index) { return *qubits[index]; }

        const std::string& getName() const { return circuit_name; }

        // --- Gate Operations ---

        void I(qubit& q) {
            validate_qubit_in_circuit(q);
            q.addGate(gate::I);
        }

        void H(qubit& q) {
            validate_qubit_in_circuit(q);
            math_H(q.index);
            q.addGate(gate::H);
        }

        void X(qubit& q) {
            validate_qubit_in_circuit(q);
            math_X(q.index);
            q.addGate(gate::X);
        }

        void Z(qubit& q) {
            validate_qubit_in_circuit(q);
            math_Z(q.index);
            q.addGate(gate::Z);
        }

        // Rk gate: applies e^(2*pi*i / 2^k) phase to the |1> component.
        // k=1 -> Z, k=2 -> S (90°), k=3 -> T (45°), etc.
        void Rk(qubit& q, int k) {
            if (k < 1)
                throw std::runtime_error("Rk: k must be >= 1");
            validate_qubit_in_circuit(q);
            math_Rk(q.index, k);
            q.addGate(gate::Rk, "", k);
        }

        void CNOT(qubit& control, qubit& target) {
            validate_qubit_in_circuit(control);
            validate_qubit_in_circuit(target);
            math_CNOT(control.index, target.index);

            size_t max_len = std::max(control.gateCount, target.gateCount);
            while (control.gateCount < max_len) I(control);
            while (target.gateCount < max_len) I(target);

            control.addGate(gate::Control);
            target.addGate(gate::Target);
        }

        void apply_circuit(const circuit& other, const std::vector<qubit*>& target_qubits) {
            if (other.size() != target_qubits.size()) {
                throw std::runtime_error("Circuit size mismatch in apply_circuit");
            }

            std::vector<qubit*> mapped_qubits(other.size());
            for (size_t i = 0; i < other.size(); ++i) {
                mapped_qubits[i] = target_qubits[i];
            }

            size_t max_gates = 0;
            for (size_t i = 0; i < other.size(); ++i) {
                max_gates = std::max(max_gates, other.qubits[i]->gateCount);
            }

            for (size_t col = 0; col < max_gates; ++col) {
                std::vector<int> controls;
                std::vector<int> targets;
                std::set<int> cnot_indices;

                for (size_t i = 0; i < other.size(); ++i) {
                    if (col < other.qubits[i]->gateCount) {
                        gate g = other.qubits[i]->getGate(col);
                        if (g == gate::Control) controls.push_back(i);
                        if (g == gate::Target)  targets.push_back(i);
                    }
                }

                std::set<int> used_targets;
                for (int ctrl_local_idx : controls) {
                    int best_target = -1;
                    int min_dist = INT_MAX;
                    for (int tgt_local_idx : targets) {
                        if (used_targets.count(tgt_local_idx)) continue;
                        int dist = std::abs(tgt_local_idx - ctrl_local_idx);
                        if (dist < min_dist) {
                            min_dist = dist;
                            best_target = tgt_local_idx;
                        }
                    }
                    if (best_target != -1) {
                        used_targets.insert(best_target);
                        cnot_indices.insert(ctrl_local_idx);
                        cnot_indices.insert(best_target);
                        math_CNOT(
                            mapped_qubits[ctrl_local_idx]->index,
                            mapped_qubits[best_target]->index
                        );
                    }
                }

                for (size_t i = 0; i < other.size(); ++i) {
                    if (cnot_indices.count(i)) continue;
                    if (col >= other.qubits[i]->gateCount) continue;

                    gate g = other.qubits[i]->getGate(col);
                    size_t idx = mapped_qubits[i]->index;

                    switch (g) {
                        case gate::H:  math_H(idx);  break;
                        case gate::X:  math_X(idx);  break;
                        case gate::Z:  math_Z(idx);  break;
                        case gate::Rk: math_Rk(idx, other.qubits[i]->getRkValue(col)); break;
                        default: break;
                    }
                }
            }

            std::string display_name =
                other.circuit_name.empty() ? " " : other.circuit_name;

            size_t global_align = 0;
            for (auto* q : qubits)
                global_align = std::max(global_align, q->gateCount);
            for (auto* q : qubits)
                while (q->gateCount < global_align)
                    q->addGate(gate::I);

            for (auto* q : target_qubits)
                q->addGate(gate::SubCircuit, display_name);

            for (auto* q : qubits)
                if (std::find(target_qubits.begin(), target_qubits.end(), q) == target_qubits.end())
                    q->addGate(gate::I);
        }

        std::string print_circuit() const {
            if (qubits.empty()) return "";

            size_t num_qubits = qubits.size();
            size_t num_rows   = num_qubits * 2 - 1;
            size_t max_gates  = 0;
            size_t max_name_len = 0;

            for (auto* q : qubits) {
                max_gates = std::max(max_gates, q->gateCount);
                std::string qname = q->name.empty() ? "q" + std::to_string(q->index) : q->name;
                max_name_len = std::max(max_name_len, qname.size());
            }

            std::vector<std::string> lines(num_rows, "");
            for (size_t i = 0; i < num_qubits; ++i) {
                std::string qname = qubits[i]->name.empty()
                    ? "q" + std::to_string(qubits[i]->index) : qubits[i]->name;
                std::string padding(max_name_len - qname.size(), ' ');
                lines[2*i] = std::to_string(qubits[i]->index) + ": " + padding + qname + " ";
                if (i < num_qubits - 1) lines[2*i + 1] = std::string(lines[2*i].size(), ' ');
            }

            for (size_t col = 0; col < max_gates; ++col) {
                // --- Subcircuit spans ---
                std::vector<int> sub_indices;
                for (size_t i = 0; i < num_qubits; ++i) {
                    auto it = qubits[i]->getSubCircuitNames().find(col);
                    if (it != qubits[i]->getSubCircuitNames().end()) sub_indices.push_back((int)i);
                }

                struct SubSpan { int start; int end; std::string name; };
                std::vector<SubSpan> subspans;
                for (size_t k = 0; k < sub_indices.size(); ) {
                    int run_start = sub_indices[k];
                    int run_end = run_start;
                    size_t j = k + 1;
                    while (j < sub_indices.size() && sub_indices[j] == (int)run_end + 1) {
                        run_end = sub_indices[j++];
                    }
                    std::string name = qubits[run_start]->getSubCircuitNames().at(col);
                    subspans.push_back({run_start, run_end, name});
                    k = j;
                }

                std::vector<bool> in_any_subspan(num_qubits, false);
                for (const auto& span : subspans)
                    for (int i = span.start; i <= span.end; ++i) in_any_subspan[i] = true;

                // --- CNOTs ---
                std::vector<int> cnot_indices;
                for (size_t i = 0; i < num_qubits; ++i) {
                    if (col < qubits[i]->gateCount) {
                        gate g = qubits[i]->getGate(col);
                        if (g == gate::Control || g == gate::Target) cnot_indices.push_back((int)i);
                    }
                }
                std::vector<std::pair<int,int>> cnot_pairs;
                for (size_t k = 0; k + 1 < cnot_indices.size(); k += 2) {
                    int a = cnot_indices[k], b = cnot_indices[k+1];
                    if (a > b) std::swap(a, b);
                    cnot_pairs.emplace_back(a, b);
                }

                // --- Column width: Rk label is "Rk(n)" so up to 6 chars + 2 dashes ---
                size_t col_width = 3;
                for (const auto& s : subspans)
                    col_width = std::max(col_width, s.name.length() + 4);
                // Check if any qubit in this column has an Rk gate — widen for label
                for (size_t i = 0; i < num_qubits; ++i) {
                    if (col < qubits[i]->gateCount && qubits[i]->getGate(col) == gate::Rk) {
                        int k = qubits[i]->getRkValue(col);
                        // label: Rk(k) e.g. "Rk(2)" = 5 chars + 2 dashes = 7 total
                        std::string label = "Rk(" + std::to_string(k) + ")";
                        col_width = std::max(col_width, label.size() + 2);
                    }
                }

                // --- Render rows ---
                for (size_t row = 0; row < num_rows; ++row) {
                    bool is_qubit_line = (row % 2 == 0);
                    size_t qubit_idx  = row / 2;

                    // Subcircuit
                    bool handled_by_subspan = false;
                    for (const auto& span : subspans) {
                        if (qubit_idx >= (size_t)span.start && qubit_idx <= (size_t)span.end) {
                            handled_by_subspan = true;
                            int canvas_sub_start = span.start * 2;
                            int canvas_sub_end   = span.end   * 2;
                            int canvas_center    = (canvas_sub_start + canvas_sub_end) / 2;

                            if (col == 0) lines[row] += is_qubit_line ? "-" : " ";

                            if (row == (size_t)canvas_sub_start) {
                                std::string border(col_width - 2, '-');
                                lines[row] += "+" + border + "+";
                            } else if (row == (size_t)canvas_sub_end && canvas_sub_start != canvas_sub_end) {
                                std::string border(col_width - 2, '-');
                                lines[row] += "+" + border + "+";
                            } else if (row == (size_t)canvas_center) {
                                size_t padding   = col_width - 2 - span.name.length();
                                size_t left_pad  = padding / 2;
                                size_t right_pad = padding - left_pad;
                                lines[row] += "|" + std::string(left_pad, ' ') + span.name
                                            + std::string(right_pad, ' ') + "|";
                            } else {
                                lines[row] += "|" + std::string(col_width - 2, ' ') + "|";
                            }
                            break;
                        }
                    }
                    if (handled_by_subspan) continue;

                    // Standard gate / spacer
                    std::string segment;
                    if (is_qubit_line) {
                        gate g = (col < qubits[qubit_idx]->gateCount)
                                 ? qubits[qubit_idx]->getGate(col) : gate::I;

                        bool is_crossing_cnot = false;
                        for (const auto& p : cnot_pairs)
                            if ((int)qubit_idx > p.first && (int)qubit_idx < p.second)
                                { is_crossing_cnot = true; break; }

                        if      (g == gate::H)       segment = "-H";
                        else if (g == gate::X)       segment = "-X";
                        else if (g == gate::Z)       segment = "-Z";
                        else if (g == gate::Control) segment = "-*";
                        else if (g == gate::Target)  segment = "-O";
                        else if (g == gate::Rk) {
                            int k = qubits[qubit_idx]->getRkValue(col);
                            segment = "-Rk(" + std::to_string(k) + ")";
                        }
                        else if (g == gate::I && is_crossing_cnot) segment = "-|";
                        else segment = "--";

                        while (segment.length() < col_width) segment += "-";
                    } else {
                        size_t upper_q = qubit_idx;
                        size_t lower_q = qubit_idx + 1;
                        bool in_cnot_path = false;
                        for (const auto& p : cnot_pairs)
                            if ((int)upper_q >= p.first && (int)lower_q <= p.second)
                                { in_cnot_path = true; break; }

                        if (in_cnot_path || (in_any_subspan[upper_q] && in_any_subspan[lower_q]))
                            segment = " | ";
                        else
                            segment = "   ";

                        while (segment.length() < col_width) segment += " ";
                    }
                    lines[row] += segment;
                }
            }

            std::string result;
            for (size_t row = 0; row < num_rows; ++row) {
                if (row % 2 == 0) lines[row] += "--";
                result += lines[row] + "\n";
            }
            return result;
        }

        std::string print_amplitudes(int precision = 2) const {
            std::string result;
            size_t n = qubits.size();
            bool first = true;

            for (size_t i = 0; i < amplitudes.size(); ++i) {
                cx amp = amplitudes[i];
                if (std::abs(amp) < 1e-12) continue;

                // Build ket label
                std::string ket;
                ket.reserve(n);
                for (size_t q = 0; q < n; ++q) {
                    size_t bit = 1ULL << (n - q - 1);
                    ket.push_back((i & bit) ? '1' : '0');
                }

                double factor = std::pow(10, precision);
                double re = std::round(amp.real() * factor) / factor;
                double im = std::round(amp.imag() * factor) / factor;

                // Format: (re + im*i) skipping zero parts for readability
                std::string amp_str;
                auto fmt = [&](double v) -> std::string {
                    std::string s = std::to_string(v);
                    s.erase(s.find_last_not_of('0') + 1, std::string::npos);
                    if (s.back() == '.') s.pop_back();
                    return s;
                };

                if (std::abs(im) < 1e-12) {
                    amp_str = fmt(re);
                } else if (std::abs(re) < 1e-12) {
                    amp_str = fmt(im) + "i";
                } else {
                    amp_str = fmt(re) + (im >= 0 ? "+" : "") + fmt(im) + "i";
                }

                if (!first) result += " + ";
                result += "(" + amp_str + ")|" + ket + ">";
                first = false;
            }
            return result.empty() ? "0" : result;
        }
};