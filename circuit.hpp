#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <climits>
#include <iostream> // Added for safety, though mostly string ops used

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

enum class gate {
    None,
    I,        // identity / placeholder
    H,        // Hadamard
    X,        // Pauli-X  (bit flip)
    Z,        // Pauli-Z  (phase flip)
    Control,  // control dot for CNOT
    Target,   // target + for CNOT
    SubCircuit // marker for applied sub-circuit
};

// Forward declaration
class circuit;

class qubit {
    private:
        std::vector<gate> gates;
        // Track sub-circuit applications: map column_index -> subcircuit_name
        std::map<size_t, std::string> subCircuitNames;
        
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
        
        void addGate(gate g, const std::string& subCircuitName = "") {
            gates.push_back(g);
            if (g == gate::SubCircuit && !subCircuitName.empty()) {
                subCircuitNames[gateCount] = subCircuitName;
            }
            gateCount++;
        }
        
        gate getGate(size_t index) const {
            if (index < gates.size()) {
                return gates[index];
            }
            return gate::None;
        }
        
        const std::map<size_t, std::string>& getSubCircuitNames() const {
            return subCircuitNames;
        }
};

class circuit {
    private:
        std::vector<qubit*> qubits;
        std::string circuit_name;

        // Ensures a given qubit belongs to this circuit.
        // Throws runtime_error if not.
        void validate_qubit_in_circuit(const qubit& q) const {
            for (auto* qp : qubits) {
                if (qp == &q) return;
            }
            throw std::runtime_error(
                "Error: Qubit \"" + q.name + "\" is not part of circuit \"" + circuit_name + "\"."
            );
        }

        
        // --- Private Math Helpers (Decoupled from Visuals) ---
        
        void math_H(size_t qubit_idx) {
            size_t n = qubits.size();
            size_t bit = 1ULL << (n - qubit_idx - 1);
            std::vector<double> new_amplitudes(amplitudes.size());
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
            std::vector<double> new_amplitudes = amplitudes;
            
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

        void math_CNOT(size_t ctrl_idx, size_t tgt_idx) {
            size_t n = qubits.size();
            size_t control_bit = 1ULL << (n - ctrl_idx - 1);
            size_t target_bit  = 1ULL << (n - tgt_idx - 1);
            std::vector<double> new_amplitudes = amplitudes;
            
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
        std::vector<double> amplitudes;
        
        circuit(const std::vector<qubit*>& q, const std::string& name = "") 
            : qubits(q), circuit_name(name) {
            amplitudes.resize(1 << qubits.size());
            std::fill(amplitudes.begin(), amplitudes.end(), 0.0);
            amplitudes[0] = 1.0;
            for(size_t i = 0; i < qubits.size(); ++i) {
                qubits[i]->index = i;
            }
        }
        
        size_t size() const {
            return qubits.size();
        }
        
        qubit& operator[](size_t index) {
            return *qubits[index];
        }
        
        const std::string& getName() const {
            return circuit_name;
        }
        
        // --- Gate Operations (Visuals + Math) ---

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
        
        void CNOT(qubit& control, qubit& target) {
            validate_qubit_in_circuit(control);
            validate_qubit_in_circuit(target);
            math_CNOT(control.index, target.index);
            
            // Visual Alignment Logic
            size_t max_len = std::max(control.gateCount, target.gateCount);
            while (control.gateCount < max_len) I(control);
            while (target.gateCount < max_len) I(target);
            
            control.addGate(gate::Control);
            target.addGate(gate::Target);
        }

        // --- Apply Subcircuit ---
        
        void apply_circuit(const circuit& other, const std::vector<qubit*>& target_qubits) {
            if (other.size() != target_qubits.size()) {
                throw std::runtime_error("Circuit size mismatch in apply_circuit");
            }
            
            // Map other's local indices to target pointers
            std::vector<qubit*> mapped_qubits(other.size());
            for (size_t i = 0; i < other.size(); ++i) {
                mapped_qubits[i] = target_qubits[i];
            }
            
            // 1. REPLAY THE MATH
            // We must simulate the 'other' circuit column by column to ensure correct order
            // relative to the internal state, but we apply the math to *this* circuit's amplitudes.
            
            size_t max_gates = 0;
            for (size_t i = 0; i < other.size(); ++i) {
                max_gates = std::max(max_gates, other.qubits[i]->gateCount);
            }
            
            for (size_t col = 0; col < max_gates; ++col) {
                // We need to identify CNOT pairs vs Single gates in 'other'
                std::vector<int> controls;
                std::vector<int> targets;
                std::set<int> cnot_indices; 
                
                // Identify CNOT components
                for (size_t i = 0; i < other.size(); ++i) {
                    if (col < other.qubits[i]->gateCount) {
                        gate g = other.qubits[i]->getGate(col);
                        if (g == gate::Control) controls.push_back(i);
                        if (g == gate::Target) targets.push_back(i);
                    }
                }

                // Apply CNOT math
                std::set<int> used_targets;
                for (int ctrl_local_idx : controls) {
                    // Find nearest available target (simplest matching logic)
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
                        
                        // Perform Math on the actual target qubits in *this* circuit
                        math_CNOT(mapped_qubits[ctrl_local_idx]->index, mapped_qubits[best_target]->index);
                    }
                }
                
                // Apply Single Gate Math
                for (size_t i = 0; i < other.size(); ++i) {
                    if (cnot_indices.count(i)) continue; // Already handled
                    
                    if (col < other.qubits[i]->gateCount) {
                        gate g = other.qubits[i]->getGate(col);
                        size_t target_idx = mapped_qubits[i]->index;
                        
                        switch (g) {
                            case gate::H: math_H(target_idx); break;
                            case gate::X: math_X(target_idx); break;
                            case gate::Z: math_Z(target_idx); break;
                            // I, Control, Target, SubCircuit handled/ignored
                            default: break;
                        }
                    }
                }
            }
            
            // 2. APPLY THE VISUAL PLACEHOLDER
            // Now that math is updated, we add ONE 'SubCircuit' gate to the visualization
            
            std::string display_name = other.circuit_name.empty() ? " " : other.circuit_name;
            
            // Align all target qubits to the same column before starting the block
            size_t align_col = 0;
            for (auto* q : target_qubits) align_col = std::max(align_col, q->gateCount);
            
            for (auto* q : target_qubits) {
                while (q->gateCount < align_col) {
                    q->addGate(gate::I);
                }
                // Add the SubCircuit marker
                q->addGate(gate::SubCircuit, display_name);
            }
        }
        
        std::string print_circuit() const {
            if (qubits.empty()) return "";

            // 1. Determine dimensions
            size_t num_qubits = qubits.size();
            size_t num_rows = num_qubits * 2 - 1; // Qubit lines + Spacer lines
            size_t max_gates = 0;
            size_t max_name_len = 0;

            for (auto* q : qubits) {
                max_gates = std::max(max_gates, q->gateCount);
                std::string qname = q->name.empty() ? "q" + std::to_string(q->index) : q->name;
                max_name_len = std::max(max_name_len, qname.size());
            }

            // 2. Initialize the lines (canvas)
            std::vector<std::string> lines(num_rows, "");

            // 3. Draw Headers (Names)
            // We append " -" to qubit lines to ensure the wire starts immediately.
            for (size_t i = 0; i < num_qubits; ++i) {
                std::string qname = qubits[i]->name.empty() ? "q" + std::to_string(qubits[i]->index) : qubits[i]->name;
                std::string padding(max_name_len - qname.size(), ' ');
                
                // Qubit line: "Index: Name -"
                lines[2 * i] = std::to_string(qubits[i]->index) + ": " + padding + qname + " ";
                
                // Spacer line: fill with spaces matching the length of the qubit line
                if (i < num_qubits - 1) {
                    lines[2 * i + 1] = std::string(lines[2 * i].size(), ' ');
                }
            }

            // 4. Iterate through circuit columns (time steps)
            for (size_t col = 0; col < max_gates; ++col) {
                // Step 4a: Analyze the column
                int sub_start_idx = -1;
                int sub_end_idx = -1;
                std::string sub_name = "";
                bool has_subcircuit = false;

                // Check for SubCircuit info
                for (size_t i = 0; i < num_qubits; ++i) {
                    auto it = qubits[i]->getSubCircuitNames().find(col);
                    if (it != qubits[i]->getSubCircuitNames().end()) {
                        if (!has_subcircuit) {
                            sub_start_idx = i;
                            sub_name = it->second;
                            has_subcircuit = true;
                        }
                        sub_end_idx = i;
                    }
                }

                // --- NEW: Detect multiple distinct CNOT pairs in this column ---
                std::vector<int> cnot_indices;
                for (size_t i = 0; i < num_qubits; ++i) {
                    if (col < qubits[i]->gateCount) {
                        gate g = qubits[i]->getGate(col);
                        if (g == gate::Control || g == gate::Target) {
                            cnot_indices.push_back((int)i);
                        }
                    }
                }
                // Pair adjacent indices into disjoint spans: (indices[0], indices[1]), (indices[2], indices[3]), ...
                std::vector<std::pair<int,int>> cnot_pairs;
                for (size_t k = 0; k + 1 < cnot_indices.size(); k += 2) {
                    int a = cnot_indices[k];
                    int b = cnot_indices[k + 1];
                    if (a > b) std::swap(a, b);
                    cnot_pairs.emplace_back(a, b);
                }
                // If cnot_indices.size() is odd, the last one is ignored (malformed column), which is safer than connecting everything.

                // Step 4b: Calculate Column Width
                size_t col_width = 3; 
                if (has_subcircuit) {
                    col_width = sub_name.length() + 4; // "|  Name  |"
                }

                // Step 4c: Render
                int canvas_sub_start = (sub_start_idx != -1) ? sub_start_idx * 2 : -1;
                int canvas_sub_end   = (sub_end_idx != -1) ? sub_end_idx * 2 : -1;
                int canvas_center    = (canvas_sub_start + canvas_sub_end) / 2;

                for (size_t row = 0; row < num_rows; ++row) {
                    bool is_qubit_line = (row % 2 == 0);
                    size_t qubit_idx = row / 2;

                    // --- SubCircuit Box Logic ---
                    if (has_subcircuit && qubit_idx >= (size_t)sub_start_idx && qubit_idx <= (size_t)sub_end_idx) {

                        // FIX: Prepend missing dash when subcircuit is at column 0 (keeps single-character alignment from your previous change)
                        if (col == 0) {
                            if (is_qubit_line) {
                                lines[row] += "-";
                            } else {
                                lines[row] += " ";
                            }
                        }

                        if (row == (size_t)canvas_sub_start) {
                            std::string border(col_width - 2, '-');
                            lines[row] += "+" + border + "+";
                        }
                        else if (row == (size_t)canvas_sub_end && canvas_sub_start != canvas_sub_end) {
                            std::string border(col_width - 2, '-');
                            lines[row] += "+" + border + "+";
                        }
                        else if (row == (size_t)canvas_center) {
                            size_t padding = col_width - 2 - sub_name.length();
                            size_t left_pad = padding / 2;
                            size_t right_pad = padding - left_pad;
                            lines[row] += "|" + std::string(left_pad, ' ') + sub_name + std::string(right_pad, ' ') + "|";
                        }
                        else if (row <= (size_t)sub_end_idx + (sub_end_idx - sub_start_idx)) {
                            lines[row] += "|" + std::string(col_width - 2, ' ') + "|";
                        }

                        continue;
                    }

                    // --- Standard Gate / Spacer Logic ---
                    std::string segment = "";
                    
                    if (is_qubit_line) {
                        gate g = (col < qubits[qubit_idx]->gateCount) ? qubits[qubit_idx]->getGate(col) : gate::I;

                        // Determine if this qubit is inside any CNOT span (exclusive); if so, crossing symbol used for I
                        bool is_crossing_cnot = false;
                        for (const auto &p : cnot_pairs) {
                            if ((int)qubit_idx > p.first && (int)qubit_idx < p.second) {
                                is_crossing_cnot = true;
                                break;
                            }
                        }

                        std::string symbol = "--";
                        if (g == gate::H) symbol = "-H";
                        else if (g == gate::X) symbol = "-X";
                        else if (g == gate::Z) symbol = "-Z";
                        else if (g == gate::Control) symbol = "-*";
                        else if (g == gate::Target) symbol = "-O";
                        else if (g == gate::I && is_crossing_cnot) symbol = "-|"; 
                        else symbol = "--";

                        segment = symbol;

                        size_t current_visual_len = segment.length();

                        // Pad with dashes based on visual length
                        while (current_visual_len < col_width) {
                            segment += "-";
                            current_visual_len++;
                        }
                    } 
                    else {
                        // Spacer line
                        size_t upper_q = qubit_idx;
                        size_t lower_q = qubit_idx + 1;
                        bool in_cnot_path = false;

                        // Only draw a vertical pipe if this spacer lies within any paired CNOT span
                        for (const auto &p : cnot_pairs) {
                            if ((int)upper_q >= p.first && (int)lower_q <= p.second) {
                                in_cnot_path = true;
                                break;
                            }
                        }

                        if (in_cnot_path) segment = " | ";
                        else segment = "   ";
                        
                        // Spacer segments are ASCII only, so length() == visual length
                        while (segment.length() < col_width) segment += " ";
                    }
                    lines[row] += segment;
                }
            }

            // 5. Finalize
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
                double amp = amplitudes[i];
                if (std::abs(amp) < 1e-12) continue;
                
                std::string ket;
                ket.reserve(n);
                for (size_t q = 0; q < n; ++q) {
                    size_t bit = 1ULL << (n - q - 1);
                    ket.push_back((i & bit) ? '1' : '0');
                }
                
                double factor = std::pow(10, precision);
                double rounded = std::round(amp * factor) / factor;
                std::string amp_str = std::to_string(rounded);
                amp_str.erase(amp_str.find_last_not_of('0') + 1, std::string::npos);
                if (amp_str.back() == '.') amp_str.pop_back();
                
                if (!first) result += " + ";
                result += "(" + amp_str + ")|" + ket + ">";
                first = false;
            }
            return result.empty() ? "0" : result;
        }
};